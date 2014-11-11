#' @name LoadMotifLibrary
#' @title Load position weight matrices.
#' @description Load the file for position weight matrices for motifs.
#' @param filename A file containing MEME format: http://meme.nbcr.net/meme/doc/meme-format.html.
#' @details This function reads the MEME formatted file containing motif information and convert them into a list of position weight matrices.
#' @return A list object of two components:
#' \tabular{ll}{
#' matrix \tab A list of position weight matrices.\cr
#' prior \tab A vector of the prior distribution for ACGT.\cr}
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{
#' pwms <- LoadMotifLibrary("/p/keles/ENCODE-CHARGE/volume1/ENCODE-Motifs/encode_motifs_for_fimo.txt")
#' }
#' @useDynLib MotifAnalysis
#' @export
LoadMotifLibrary <- function(filename) {
    lines <- readLines(filename)
    bkngLineNum <- grep("Background letter frequencies", lines) + 1
    priorProb <- as.numeric(strsplit(lines[bkngLineNum], " ")[[1]][seq(4) * 2])
    motifLineNums <- grep("MOTIF", lines)
    motifnames <-
        sapply(strsplit(lines[motifLineNums], " "), function(x) x[2])
    allmats <- as.list(seq_along(motifnames))

    for(matrixId in seq_along(motifLineNums)) {
        motifLineNum <- motifLineNums[matrixId]
        nrows <- as.integer(strsplit(lines[motifLineNum + 1], " ")[[1]][6])
        pwm <-
          t(matrix(as.numeric(unlist(strsplit(lines[seq(nrows) + motifLineNum + 1], " "))), nrow = 4))
        pwm <- t(apply(pwm, 1,
                     function(x) {
                     x[x < 1e-10] <- 1e-10 / (1 - 1e-10 * sum(x < 1e-10)) * sum(x)
                     return(x / sum(x))
                   }))
        allmats[[matrixId]] <- pwm
    }
    names(allmats) <- motifnames
    return(
        list(matrix = allmats,
             prior = priorProb))
}


#' @name LoadSNPData
#' @title Load the SNP information and code the genome sequences around the SNP locations.
#' @description Load the SNP data.
#' @param filename A table containing the SNP information. Must contain at least four columns:
#' \tabular{ll}{
#' chr \tab chromosome.\cr
#' snp \tab The nucleotide position of the SNP.\cr
#' snpid \tab The names of the SNPs.\cr
#' a1 \tab The deoxyribose at the reference genome.\cr
#' a2 \tab The deoxyribose at the SNP genome.\cr
#' }
#' @param genome.lib A string of the library name for the genome version. Default: "BSgenome.Hsapiens.UCSC.hg19".
#' @param half.window.size An integer for the half window size around the SNP within which the motifs are matched. Default: 30.
#' @param ... Other parameters passed to 'read.table'.
#' @details TODO.
#' @return A list object containing the following components:
#' \tabular{ll}{
#' sequence_matrix \tab A list of integer vectors representing the deroxyribose sequence around each SNP.\cr
#' a1 \tab An integer vector for the deroxyribose at the SNP location on the reference genome.\cr
#' a2 \tab An integer vector for the deroxyribose at the SNP location on the SNP genome.\cr
#' }
#' The results are coded as: "A"-1, "C"-2, "G"-3, "T"-4.
#' @author Chandler Zuo \email{zuo@@stat.wisc.edu}
#' @examples
#' \dontrun{LoadSNPData("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed")}
#' @useDynLib MotifAnalysis
#' @export
LoadSNPData <- function(filename, genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
                        half.window.size = 30, ...) {
  ## load the corresponding genome version
  library(package = genome.lib, character.only = TRUE)
  tbl <- read.table(filename, header = TRUE, stringsAsFactors = FALSE, ...)
  ## check if the input file has the required information
  if(sum(!c("snp", "chr", "a1", "a2", "snpid") %in% names(tbl)) > 0) {
    stop("Error: 'filename' must be a table containing 'snp' and 'chr' columns.")
  }
  seqvec <- getSeq(Hsapiens,
                   as.character(tbl$chr),
                   start=tbl$snp - half.window.size,
                   end=tbl$snp + half.window.size,
                   as.character=TRUE)
  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  sequences <- sapply(seqvec, function(x) codes[strsplit(x, "")[[1]]])
  colnames(sequences) <- tbl$snpid
  rownames(sequences) <- NULL
  a1 <- codes[tbl$a1]
  a2 <- codes[tbl$a2]
  names(a1) <- names(a2) <- NULL
  transition <- .Call("transition_matrix", sequences, package = "MotifAnalysis")
  prior <- apply(transition, 1, sum)
  prior <- prior / sum(prior)
  transition <- transition / apply(transition, 1, sum)
  names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  a1.ref.base.id <- which(a1 == sequences[half.window.size + 1, ])
  a2.ref.base.id <- which(a2 == sequences[half.window.size + 1, ])
  ## store SNPs that have the same base in SNP and REF alleles only once
  a2.ref.base.id <- a2.ref.base.id[!a2.ref.base.id %in% a1.ref.base.id]
  sequences <- sequences[, c(a1.ref.base.id, a2.ref.base.id)]
  ref.base <- c(a1[a1.ref.base.id], a2[a2.ref.base.id])
  snp.base <- c(a2[a1.ref.base.id], a1[a2.ref.base.id])
  return(list(
              sequence_matrix= sequences,
              ref_base = ref.base,
              snp_base = snp.base,
              transition = transition,
              prior = prior
              ))
}

#' @name ComputeMotifScore
#' @title Compute the scores for SNP effects on motifs.
#' @description Compute the log-likelihood scores for motifs.
#' @param motif.lib A list object with the output format of function 'LoadMotifLibrary'.
#' @param snp.info A list object with the output format of function 'LoadSNPData'.
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @details TODO.
#' @return A data.table containing:
#' \tabular{ll}{
#' motif \tab Name of the motif.\cr
#' snpid \tab The SNP id.\cr
#' ref_seq \tab The nucleobase sequence for the reference allele.\cr
#' snp_seq \tab The nucleobase sequence for the SNP allele.\cr
#' ref_seq_rev \tab The nucleobase sequence for the reference allele on the reverse strand.\cr
#' snp_seq_rev \tab The nucleobase sequence for the SNP allele on the reverse strand.\cr
#' motif_len \tab Length of the motif.\cr
#' ref_start, ref_end, ref_strand \tab Location of the best matching subsequence on the reference allele.\cr
#' snp_start, snp_end, snp_strand \tab Location of the best matching subsequence on the SNP allele.\cr
#' log_lik_ref \tab Log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab Log-likelihood score for the SNP allele.\cr
#' log_lik_ratio \tab The log-likelihood ratio.\cr
#' log_enhance_odds \tab Difference in log-likelihood ratio between SNP allele and reference allele based on the best matching subsequence on the reference allele.\cr
#' log_reduce_odds \tab Difference in log-likelihood ratio between reference allele and SNP allele based on the best matching subsequence on the SNP allele.\cr
#' ref_match_seq \tab Best matching subsequence on the reference allele.\cr
#' snp_match_seq \tab Best matching subsequence on the SNP allele.\cr
#' ref_match_seq \tab Subsequence on the reference allele corresponding to the best matching location on the SNP allele.\cr
#' snp_match_seq \tab Subsequence on the SNP allele corresponding to the best matching location on the reference allele.\cr
#' }
#' Each component is a matrix, with columns representing the motifs and the rows representing the SNPs.
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' ComputeMotifScore(motif_library, snpInfo, ncores = 2)
#' @useDynLib MotifAnalysis
#' @import data.table doMC
#' @export
ComputeMotifScore <- function(motif.lib, snp.info, ncores = 1) {
  ## check arguments
  if(sum(!c("prior", "matrix") %in% names(motif.lib)) > 0) {
    stop("Error: 'motif.lib' must contain both 'matrix' and'prior'.")
  }
  if(sum(!unlist(sapply(motif.lib$matrix, is.matrix))) > 0 | sum(unlist(sapply(motif.lib$matrix, ncol)) != 4) > 0) {
    stop("Error: 'motif.lib$matrix' must be a list of numeric matrices each with 4 columns.")
  }
  if(sum(!c("sequence_matrix", "snp_base", "ref_base") %in% names(snp.info)) > 0) {
    stop("Error: 'snp.info' must contain three components: 'ref_base', 'snp_base', 'sequence_matrix'.")
  }
  if(ncol(snp.info$sequence_matrix) != length(snp.info$ref_base) | length(snp.info$ref_base) != length(snp.info$snp_base)) {
    stop("Error: the number of columns of 'snp.info$sequence_matrix', the length of 'snp.info$ref_base' and the length of 'snp.info$snp_base' must be the same.")
  }
  if(sum(sort(unique(c(c(snp.info$sequence_matrix), snp.info$ref_base, snp.info$snp_base))) != seq(4)) > 0) {
    stop("Error: 'snp.info$sequence_matrix', 'snp.info$ref_base', 'snp.info$snp_base' can only contain entries in 1, 2, 3, 4.")
  }
  if(nrow(snp.info$sequence_matrix) / 2 == as.integer(nrow(snp.info$sequence_matrix) / 2)) {
    stop("Error: 'snp.info$sequence_matrix' must have an odd number of rows so that the central row refers to the SNP nucleotide.")
  }
  
  registerDoMC(ncores)

  motif_score_par <- foreach(i = seq(ncores)) %dopar% {
    k <- as.integer(length(snp.info$ref_base) / ncores)
    if(i < ncores) {
      ids <- seq(k) + k * (i - 1)
    } else {
      ids <- (k * (ncores - 1) + 1):length(snp.info$ref_base) 
    }
    this.snp.info <- list(sequence_matrix = snp.info$sequence_matrix[, ids],
                          ref_base = snp.info$ref_base[ids], snp_base = snp.info$snp_base[ids])
    .Call("motif_score", motif.lib, this.snp.info, package = "MotifAnalysis")
  }

  motif.scores <- motif_score_par[[1]]
  if(ncores > 1) {
    for(i in 2:ncores) {
      for(j in seq_along(motif.scores)) {
        motif.scores[[j]] <- rbind(motif.scores[[j]], motif_score_par[[i]][[j]])
      }
    }
  }

  motifs <- names(motif.lib$matrix)
  snpids <- colnames(snp.info$sequence_matrix)
  nsnps <- ncol(snp.info$sequence_matrix)
  nmotifs <- length(motif.lib$matrix)
  for(i in seq_along(motif.scores)) {
    rownames(motif.scores[[i]]) <- snpids
    colnames(motif.scores[[i]]) <- motifs
  }


  motif_tbl <- data.table(motif = motifs,
                          motif_len = sapply(motif.lib$matrix, nrow))

  ## sequences on the reference genome
  ref_seqs <- apply(snp.info$sequence_matrix, 2, function(x) paste(c("A", "C", "G", "T")[x], collapse = ""))
  ref_seqs_rev <- apply(snp.info$sequence_matrix, 2, function(x) paste(c("A", "C", "G", "T")[5 - rev(x)], collapse = ""))
  ## sequences on the snp allele
  id1 <- seq(as.integer(nrow(snp.info$sequence_matrix) / 2))
  id2 <- id1 + (nrow(snp.info$sequence_matrix) + 1) / 2
  snp_seqs <- apply(rbind(snp.info$sequence_matrix, snp.info$snp_base), 2,
                    function(x)
                    paste(c("A", "C", "G", "T")[x[c(id1, length(x), id2)]],
                          collapse = "")
                    )
  snp_seqs_rev <- apply(rbind(snp.info$sequence_matrix, snp.info$snp_base), 2,
                    function(x)
                    paste(c("A", "C", "G", "T")[5 - rev(x[c(id1, length(x), id2)])],
                          collapse = "")
                    )

  snp_tbl <- data.table(snpid = snpids,
                        ref_seq = ref_seqs,
                        snp_seq = snp_seqs,
                        ref_seq_rev = ref_seqs_rev,
                        snp_seq_rev = snp_seqs_rev)

  len_seq <- nrow(snp.info$sequence_matrix)
  strand_ref <- (motif.scores$match_ref_base > 0)
  ref_start <- motif.scores$match_ref_base
  ref_start[!strand_ref] <- len_seq + motif.scores$match_ref_base[!strand_ref] + 1

  strand_snp <- (motif.scores$match_snp_base > 0)
  snp_start <- motif.scores$match_snp_base
  snp_start[!strand_snp] <- len_seq + motif.scores$match_snp_base[!strand_snp] + 1

  motif_score_tbl <- data.table(snpid = rep(snpids, nmotifs),
                                 motif = rep(motifs, each = nsnps),
                                 log_lik_ref = c(motif.scores$log_lik_ref),
                                 log_lik_snp = c(motif.scores$log_lik_snp),
                                 log_lik_ratio = c(motif.scores$log_lik_ratio),
                                 log_enhance_odds = c(motif.scores$log_enhance_odds),
                                 log_reduce_odds = c(motif.scores$log_reduce_odds),
                                 ref_start = c(ref_start),
                                 snp_start = c(snp_start),
                                 ref_strand = c("-", "+")[1 + as.integer(strand_ref)],
                                 snp_strand = c("-", "+")[1 + as.integer(strand_snp)]
                                 )

  setkey(motif_score_tbl, motif)
  setkey(motif_tbl, motif)
  motif_score_tbl <- motif_tbl[motif_score_tbl]
  motif_score_tbl[ref_strand == "-", ref_start := ref_start - motif_len + 1]
  motif_score_tbl[, ref_end := ref_start + motif_len - 1]
  motif_score_tbl[snp_strand == "-", snp_start := snp_start - motif_len + 1]
  motif_score_tbl[, snp_end := snp_start + motif_len - 1]

  setkey(motif_score_tbl, snpid)
  setkey(snp_tbl, snpid)
  motif_score_tbl <- snp_tbl[motif_score_tbl]
  motif_score_tbl[ref_strand == "+", ref_match_seq := substr(ref_seq, ref_start, ref_end)]
  motif_score_tbl[ref_strand == "-", ref_match_seq := substr(ref_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
  motif_score_tbl[snp_strand == "+", snp_match_seq := substr(snp_seq, snp_start, snp_end)]
  motif_score_tbl[snp_strand == "-", snp_match_seq := substr(snp_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]
  motif_score_tbl[ref_strand == "+", snp_seq_ref_match := substr(snp_seq, ref_start, ref_end)]
  motif_score_tbl[ref_strand == "-", snp_seq_ref_match := substr(snp_seq_rev, len_seq - ref_end + 1, len_seq - ref_start + 1)]
  motif_score_tbl[snp_strand == "+", ref_seq_snp_match := substr(ref_seq, snp_start, snp_end)]
  motif_score_tbl[snp_strand == "-", ref_seq_snp_match := substr(ref_seq_rev, len_seq - snp_end + 1, len_seq - snp_start + 1)]
  
  return(motif_score_tbl[, list(snpid,
                                motif,
                                ref_seq,
                                snp_seq,
                                ref_seq_rev,
                                snp_seq_rev,
                                motif_len,
                                ref_start,
                                ref_end,
                                ref_strand,
                                snp_start,
                                snp_end,
                                snp_strand,
                                log_lik_ref,
                                log_lik_snp,
                                log_lik_ratio,
                                log_enhance_odds,
                                log_reduce_odds,
                                ref_match_seq,
                                snp_match_seq,
                                ref_seq_snp_match,
                                snp_seq_ref_match)])
  
}

CheckSameLength <- function(x) {
  if(length(x) == 1) {
    return(TRUE)
  }
  return(var(unlist(sapply(x, length))) == 0)
}

#' @name ComputePValues
#' @title Compute p values.
#' @description TODO.
#' @param motif.lib A list object with the output format of function 'LoadMotifLibrary'.
#' @param snp.info A list object with the output format of function 'LoadSNPData'.
#' @param motif.scores A data.table object containing at least the following columns:
#' \tabular{ll}{
#' motif \tab The name of the motif.\cr
#' log_lik_ref \tab The log-likelihood score for the reference allele.\cr
#' log_lik_snp \tab The log-likelihood score for the SNP allele.\cr
#' }
#' @param ncores An integer for the number of parallel process. Default: 1.
#' @details TODO.
#' @return A data.table extending 'motif.scores' by the following additional columns:
#' \tabular{ll}{
#' pval_ref \tab P values for scores on the reference allele.\cr
#' pval_snp \tab P values for scores on the SNP allele.\cr
#' pval_diff \tab P values for the difference in scores between the reference and the SNP alleles.\cr
#' }
#' @author Chandler Zuo\email{zuo@@stat.wisc.edu}
#' @examples
#' data(example)
#' ComputePValues(motif_library, snpInfo, motif_scores, ncores = 4) 
#' @import doMC Rcpp data.table
#' @useDynLib MotifAnalysis
#' @export
ComputePValues <- function(motif.lib, snp.info, motif.scores, ncores = 1) {
  registerDoMC(ncores)
  results <- foreach(i = seq_along(motif.lib$matrix)) %dopar% {
      rowids <- which(motif.scores$motif == names(motif.lib$matrix)[i])
    scores <- cbind(motif.scores$log_lik_ref[rowids],
                    motif.scores$log_lik_snp[rowids])
    pwm <- motif.lib$matrix[[i]]
    pwm[pwm < 1e-10] <- 1e-10
    wei.mat <- pwm
    for(i in seq(nrow(wei.mat))) {
        for(j in seq(ncol(wei.mat))) {
            wei.mat[i, j] <- exp(mean(log(pwm[i, j] / pwm[i, -j])))
        }
    }
    pval_a <- .Call("test_p_value_diff", pwm, snp.info$prior,
                    snp.info$transition, scores, 0.01,
                    package = "MotifAnalysis")
    pval_diff_r <- .Call("test_p_value_diff", pwm,
                         wei.mat, pwm ^ 0.5, snp.info$prior,
                         snp.info$transition, scores,
                         0.1, package = "MotifAnalysis")
    message("Finished testing the ", i, "th motif")
      motif.scores[rowids, pval_ref := pval_a[, 1]]
      motif.scores[rowids, pval_snp := pval_a[, 2]]
      motif.scores[rowids, pval_diff := pval_diff_r]
  }
  return(motif.scores)
}
