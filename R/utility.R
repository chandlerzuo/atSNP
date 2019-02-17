.structure <- function(pval_mat) {
    if(is.matrix(pval_mat)) {
      if(nrow(pval_mat) > 1) {
        id <- apply(pval_mat[, c(2, 4)], 1, which.min)
        return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
                     pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
      }
    }
    
    pval_mat <- as.vector(pval_mat)
    id <- which.min(pval_mat[c(2, 4)])
    return(cbind(pval_mat[c(1, 3)][id],
                 pval_mat[c(2, 4)][id]))
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}

myStrSplit <- function(x, split) {
  ret <- list(seq_along(x))
  for(i in seq_along(x)) {
    ret[[i]] <- x[i]
    for(sp in split) {
      ret[[i]] <- unlist(strsplit(ret[[i]], sp))
      ret[[i]] <- ret[[i]][nchar(ret[[i]]) > 0]
      if(length(ret[[i]]) == 0)
        break
    }
  }
  return(ret)
}

checkSNPids<-function(ids) {
  return(all(any(!is(ids, "character"),  length(is)==0), is.null(is)==FALSE))
  
}
checkMotifs<-function(m) {
  return(all(any(!is(m, "character"),  length(m)==0), is.null(m)==FALSE))
}

motif_score_par<-function(i, par.k=k, par.ncores=ncores, par.motifs=motifs, par.nmotifs=nmotifs, par.snpids=snpids,  par.snpbases=snpbases, par.len_seq=len_seq, par.motif.lib=motif.lib, par.snp.info=snp.info) {
  if(i < par.ncores) {
    ids <- seq(par.k) + par.k * (i - 1)
  } else {
    ids <- (par.k * (par.ncores - 1) + 1):length(par.snp.info$ref_base)
  }
  this.snp.info <- list(sequence_matrix = t(t(par.snp.info$sequence_matrix[, ids])),
                        ref_base = par.snp.info$ref_base[ids], snp_base = par.snp.info$snp_base[ids])
  motif.scores <- .Call("motif_score", par.motif.lib, this.snp.info, package = "atSNP")
  for(j in seq_along(motif.scores)) {
    rownames(motif.scores[[j]]) <- par.snpids[ids]
    colnames(motif.scores[[j]]) <- par.motifs
  }
  
  nids<-length(ids)
  motif_tbl <- data.table(motif = par.motifs, motif_len = sapply(par.motif.lib, nrow))
  motif_len.m<-matrix(motif_tbl$motif_len, nids, par.nmotifs, byrow=TRUE)
  rownames(motif_len.m)<-par.snpids[ids]
  
  strand_ref <- (motif.scores$match_ref_base > 0)
  ref_start <- motif.scores$match_ref_base
  ref_start[!strand_ref] <- par.len_seq + motif.scores$match_ref_base[!strand_ref] - motif_len.m[!strand_ref]+2
  ref_end <- ref_start + motif_len.m - 1
  
  strand_snp <- (motif.scores$match_snp_base > 0)
  snp_start <- motif.scores$match_snp_base
  snp_start[!strand_snp] <- par.len_seq + motif.scores$match_snp_base[!strand_snp] - motif_len.m[!strand_snp]+2
  snp_end <- snp_start + motif_len.m - 1
  motif_score_tbl_dt <- data.table(snpid = rep(par.snpids[ids], par.nmotifs),
                                   motif = rep(par.motifs, each = length(ids)),
                                   log_lik_ref = c(motif.scores$log_lik_ref),
                                   log_lik_snp = c(motif.scores$log_lik_snp),
                                   log_lik_ratio = c(motif.scores$log_lik_ratio),
                                   log_enhance_odds = c(motif.scores$log_enhance_odds),
                                   log_reduce_odds = c(motif.scores$log_reduce_odds),
                                   ref_start = c(ref_start),
                                   snp_start = c(snp_start),
                                   ref_end=c(ref_end),
                                   snp_end=c(snp_end),
                                   ref_strand = c("-", "+")[1 + as.integer(strand_ref)],
                                   snp_strand = c("-", "+")[1 + as.integer(strand_snp)],
                                   snpbase= rep(par.snpbases[ids], par.nmotifs)                                                             
  )
  setkey(motif_score_tbl_dt, motif, snpbase)
  setkey(motif_tbl, motif)
  motif_score_tbl <- as.data.frame(motif_tbl[motif_score_tbl_dt], stringAsFactors = FALSE)
  return(motif_score_tbl)
  
}

match_subseq_par<-function(i, par.k=k, par.ncores=ncores, par.snp.tbl=snp.tbl, par.snpids=snpids, par.motif.scores=motif.scores, par.motif=motif, par.motif.tbl=motif.tbl) {
  if(i < par.ncores) {
    ids <- seq(par.k) + par.k * (i - 1)
  } else {
    ids <- (par.k * (par.ncores - 1) + 1):length(par.snpids)
  }
  motif.scores_i <- par.motif.scores[snpid %in% par.snpids[ids], ]
  setkey(motif.scores_i, motif)
  motif.scores_i <- par.motif.tbl[motif.scores_i]
  setkey(motif.scores_i, snpid, snpbase)
  motif.scores_i <- par.snp.tbl[motif.scores_i]
  motif.scores_i[ref_strand == "+", ref_match_seq := substr(ref_seq, ref_start, ref_end)]
  motif.scores_i[ref_strand == "-", ref_match_seq := substr(ref_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
  motif.scores_i[snp_strand == "+", snp_match_seq := substr(snp_seq, snp_start, snp_end)]
  motif.scores_i[snp_strand == "-", snp_match_seq := substr(snp_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]
  motif.scores_i[ref_strand == "+", snp_seq_ref_match := substr(snp_seq, ref_start, ref_end)]
  motif.scores_i[ref_strand == "-", snp_seq_ref_match := substr(snp_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
  motif.scores_i[snp_strand == "+", ref_seq_snp_match := substr(ref_seq, snp_start, snp_end)]
  motif.scores_i[snp_strand == "-", ref_seq_snp_match := substr(ref_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]
  
  return(motif.scores_i[, list(snpid, motif, ref_seq, snp_seq, motif_len, ref_start, ref_end, ref_strand, snp_start,
                               snp_end, snp_strand,log_lik_ref, log_lik_snp, log_lik_ratio, log_enhance_odds, log_reduce_odds,
                               IUPAC, ref_match_seq, snp_match_seq, ref_seq_snp_match, snp_seq_ref_match, snpbase)])
}

results_motif_par<-function(i, par.prior=snp.info$prior, par.transition=snp.info$transition, par.motif.lib=motif.lib, par.motif.scores=motif.scores, par.testing.mc=testing.mc, par.figdir=figdir) {
  rowids <- which(par.motif.scores$motif == names(par.motif.lib)[i])
  scores <- cbind(par.motif.scores$log_lik_ref[rowids],
                  par.motif.scores$log_lik_snp[rowids])
  pwm <- par.motif.lib[[i]]
  pwm[pwm < 1e-10] <- 1e-10
  wei.mat <- pwm
  for(i in seq(nrow(wei.mat))) {
    for(j in seq(ncol(wei.mat))) {
      wei.mat[i, j] <- exp(mean(log(pwm[i, j] / pwm[i, -j])))
    }
  }
  
  if(nrow(scores) > 5000) {
    p <- 5 / nrow(scores)
  } else {
    p <- 1 / nrow(scores)
  }
  
  m <- 20
  b <- (1 / p) ^ ( 1 / sum(seq(m)))
  allp <- rep(1, m + 1)
  step <- b
  for(k in rev(seq(m))) {
    allp[k] <- allp[k + 1] / step
    step <- step * b
  }
  allp <- allp[-(m + 1)]
  
  score.p <- unique(quantile(c(scores), 1 - allp))
  if(length(score.p) > length(unique(c(scores)))) {
    score.p <- rev(unique(sort(c(scores))))
    allp <- seq_along(score.p) / length(c(scores))
  }
  
  pval_a <- pval_cond <- matrix(1, nrow = nrow(scores), ncol = 4)
  for(l in seq_along(allp)) {
    if(l == 1) {
      score.upp = max(scores) + 1
    } else if(l <= length(allp)) {
      score.upp = score.p[l - 1]
    } else {
      score.upp = quantile(c(scores), 0.2)
    }
    if(l >= length(allp)) {
      score.low = min(scores) - 1
    } else {
      score.low = score.p[l + 1]
    }
    compute.id <- which(scores < score.upp & scores >= score.low)
    if(length(compute.id) == 0) {
      next
    }
    if(l < length(allp) + 1) {
      theta <- .Call("test_find_theta", pwm, par.prior, par.transition, score.p[l], package = "atSNP")
    } else {
      theta <- 0
    }
    ## set the importance sample size
    if(par.testing.mc==FALSE) {
      n_sample <- 2000
      if(l <= length(allp)) {
        n_sample <- as.integer((1 - allp[l]) / allp[l] * 100)
      }
      if(n_sample > 1e5) {
        n_sample <- 1e5
      }
      if(n_sample < 5000) {
        n_sample <- 2000
      }
      pval_a.new <- .Call("test_p_value", pwm, par.prior, par.transition, scores[compute.id], theta, n_sample, package = "atSNP")
    } else {
      pval_a.new <- .Call("test_p_value", pwm, par.prior, par.transition, scores[compute.id], theta, 100, package = "atSNP")
    }
    
    pval_cond.new <- .structure(pval_a.new[, 4 + seq(4)])
    pval_a.new <- .structure(pval_a.new[, seq(4)])
    
    update.id <- which(pval_a.new[, 2] < pval_a[, 3:4][compute.id])
    pval_a[compute.id[update.id]] <- pval_a.new[update.id, 1]
    pval_a[compute.id[update.id] + 2 * nrow(pval_a)] <- pval_a.new[update.id, 2]
    update.id <- which(pval_cond.new[, 2] < pval_cond[, 3:4][compute.id])
    pval_cond[compute.id[update.id]] <- pval_cond.new[update.id, 1]
    pval_cond[compute.id[update.id] + 2 * nrow(pval_cond)] <- pval_cond.new[update.id, 2]
  }
  
  pval_a[pval_a[, seq(2)] > 1] <- 1
  pval_cond[pval_cond[, seq(2)] > 1] <- 1
  adjusted <- FALSE
  ## Force the p-values to be increasing
  while(!adjusted) {
    pval_a.sorted <- sort(pval_a[, seq(2)])[rank(-c(scores))]
    pval_cond.sorted <- sort(pval_cond[, seq(2)])[rank(-c(scores))]
    flag1 <- flag2 <- TRUE
    if(prod(pval_a.sorted == pval_a[, seq(2)]) != 1 |
       prod(pval_cond.sorted == pval_cond[, seq(2)]) != 1) {
      pval_a[, seq(2)] <- pval_a.sorted
      pval_cond[, seq(2)] <- pval_cond.sorted
      flag1 <- FALSE
    }
    ## force the conditional p-value <= p-value
    adjust.id <- which(pval_a[, seq(2)] < pval_cond[, seq(2)])
    if(length(adjust.id) > 0) {
      flag2 <- FALSE
      pval_cond[adjust.id] <- (pval_cond[adjust.id] + pval_a[adjust.id]) / 2
      pval_a[adjust.id] <- pval_cond[adjust.id]
    }
    adjusted <- flag1 & flag2
  }
  
  rank_ratio <- abs(log(pval_a[, 1] + 1e-10) - log(pval_a[, 2] + 1e-10))
  score_diff <- apply(scores, 1, function(x) abs(diff(x)))
  score.p <- round(quantile(score_diff, c((seq(8) + 1) / 10, 0.9 + seq(9) / 100)))
  if(round(quantile(score_diff, 0.1) + 1) < round(quantile(score_diff, 0.9))) {
    score.p <- c(score.p,
                 seq(round(quantile(score_diff, 0.1) + 1),
                     round(quantile(score_diff, 0.9)),
                     by = 2)
    )
  }
  score.p <- rev(sort(unique(score.p)))
  
  pval_diff <- pval_rank <- matrix(1, nrow = length(score_diff), ncol = 2)
  for(l in seq_along(score.p)) {
    if(l == 1) {
      score.upp <- max(score_diff) + 1
    } else {
      score.upp <- score.p[l - 1]
    }
    if(l == length(score.p)) {
      score.low <- min(score_diff) - 1
    } else {
      score.low <- score.p[l + 1]
    }
    compute.id <- which(score_diff < score.upp & score_diff >= score.low)
    if(length(compute.id) == 0) {
      next
    }
    ## set the importance sample size
    if(par.testing.mc==FALSE) {
      n_sample <- 2000
      p <- mean(score_diff >= score.p[l])
      if(p == 0) {
        p <- 1e-4
      }
      if(l <= length(allp)) {
        n_sample <- as.integer((1 - p) / p * 100)
      }
      if(n_sample > 1e5) {
        n_sample <- 1e5
      }
      if(n_sample < 2000) {
        n_sample <- 2000
      }
      
      pval_diff.new <- .Call("test_p_value_change", pwm,
                             wei.mat, pwm + 0.25,  par.prior,
                             par.transition, score_diff[compute.id],
                             rank_ratio[compute.id],
                             score.p[l], n_sample, package = "atSNP")
    } else {
      pval_diff.new <- .Call("test_p_value_change", pwm,
                             wei.mat, pwm + 0.25,  par.prior,
                             par.transition, score_diff[compute.id],
                             rank_ratio[compute.id],
                             score.p[l], 100, package = "atSNP")         
    }
    pval_rank.new <- .structure(pval_diff.new$rank)
    pval_diff.new <- .structure(pval_diff.new$score)
    update.id <- which(pval_diff.new[, 2] < pval_diff[compute.id, 2])
    pval_diff[compute.id[update.id], ] <- pval_diff.new[update.id, ]
    update.id <- which(pval_rank.new[, 2] < pval_rank[compute.id, 2])
    pval_rank[compute.id[update.id], ] <- pval_rank.new[update.id, ]
    ## print(summary(pval_diff.new[,1]))
  }
  
  ## force the monotonicity
  pval_diff[, 1] <- sort(pval_diff[,1])[rank(-score_diff)]
  pval_diff[pval_diff[, 1] > 1, 1] <- 1
  pval_rank[, 1] <- sort(pval_rank[,1])[rank(-rank_ratio)]
  pval_rank[pval_rank[, 1] > 1, 1] <- 1
  message("Finished testing motif No. ", i)
  if(!is.null(par.figdir)) {
    if(!file.exists(par.figdir)) {
      dir.create(par.par.figdir)
    }
    plotdat <- data.frame(
      score = c(scores),
      p.value = c(pval_a[, seq(2)]),
      var = c(pval_a[, 3:4]),
      Allele = rep(c("ref", "snp"), each = nrow(scores))
    )
    plotdat.diff <- data.frame(
      score = score_diff,
      p.value = pval_diff[,1],
      var = pval_diff[,2]
    )
    plotdat <- unique(plotdat)
    plotdat.diff <- unique(plotdat.diff)
    localenv <- environment()
    options(warn = -1)
    pdf(file.path(par.figdir, paste("motif", i, ".pdf", sep = "")), width = 10, height = 10)
    id <- which(rank(plotdat$p.value[plotdat$Allele == "ref"]) <= 500)
    print(ggplot(aes(x = score, y = p.value), data = plotdat[plotdat$Allele == "ref", ], environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(par.motif.lib)[i], "ref")))
    id <- which(rank(plotdat$p.value[plotdat$Allele == "snp"]) <= 500)
    print(ggplot(aes(x = score, y = p.value), data = plotdat[plotdat$Allele == "snp", ], environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(par.motif.lib)[i], "SNP")))
    print(ggplot(aes(x = score, y = p.value), data = plotdat.diff, environment = localenv) + geom_point() + scale_y_log10(breaks = 10 ^ seq(-8, 0)) + geom_errorbar(aes(ymax = p.value + sqrt(var), ymin = p.value - sqrt(var))) + ggtitle(paste(names(par.motif.lib)[i], " Change")))
    dev.off()
  }
  
  return(list(rowids = rowids,
              pval_a = pval_a,
              pval_cond = pval_cond,
              pval_diff = pval_diff,
              pval_rank = pval_rank))
}
