library(atSNP)
library(testthat)

if(FALSE) {
  ## Import the ENCODE motif library
  
  encode_motif <- LoadMotifLibrary("http://compbio.mit.edu/encode-motifs/motifs.txt", tag = ">", transpose = FALSE, field = 1, sep = c("\t", " ", ">"), skipcols = 1, skiprows = 1, pseudocount = 0)

  lines <- readLines("http://compbio.mit.edu/encode-motifs/motifs.txt")
  title.no <- grep(">", lines)
  source("~/atsnp_git/atSNP/R/utility.R")
  title.info <- sapply(lines[title.no], function(x) myStrSplit(x, c(">", " ", "\t")))
  nfields <- sapply(title.info, length)
  allnames <- encode_motifinfo <- rep("", length(title.no))
  for(i in 1:2) {
    allnames[nfields == i + 1] <- sapply(title.info[nfields == i + 1], function(x) x[i])
    encode_motifinfo[nfields == i + 1] <- sapply(title.info[nfields == i + 1], function(x) x[i + 1])
  }

  names(encode_motifinfo) <- names(encode_motif) <- allnames

  system.time(save(encode_motif, encode_motifinfo, file = "~/atsnp_git/atSNP/data/encode_library.rda"))
  
  ## Import the JASPAR library
  jaspar_motif <- LoadMotifLibrary(
                                   "http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt",
                                   tag = ">", skiprows = 1, skipcols = 0, transpose = TRUE, field = 1, 
                           sep = c(">", "\t", " "), pseudocount = 1)
  lines <- readLines("http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_all.txt")
  title.no <- grep(">", lines)
  source("~/atsnp_git/atSNP/R/utility.R")
  title.info <- sapply(lines[title.no], function(x) myStrSplit(x, c(">", " ", "\t")))
  nfields <- sapply(title.info, length)
  allnames <- jaspar_motifinfo <- rep("", length(title.no))
  allnames <- sapply(title.info, function(x) x[1])
  jaspar_motifinfo <- sapply(title.info, function(x) x[2])
  names(jaspar_motifinfo) <- names(jaspar_motif) <- allnames
  system.time(save(jaspar_motifinfo, jaspar_motif, file = "~/atsnp_git/atSNP/data/jaspar_library.rda"))
  
  # construct the test data set
  motif_library <- encode_motif
  system.time(snpInfo <- LoadSNPData("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed", nrow = 20))
  motif_library <- motif_library[c(1:2)]
  snp_tbl <- read.table("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed", nrow = 20, header = TRUE)[, c("snpid", "a1", "a2", "chr", "snp")]
  motif_scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 2)

  motif_pval <- ComputePValues(motif_library, snpInfo, motif_scores$motif.scores, ncores = 2)
  
  i <- 2

  motif_pval[, pval_ratio := abs(log(pval_ref + 1e-10) - log(pval_snp + 1e-10))]
  
  par(mfrow = c(2, 2))
  plot(log(pval_diff) ~ abs(log_lik_ratio), data = motif_pval[motif == names(motif_library)[i], ])
  plot(log(pval_rank) ~ pval_ratio, data = motif_pval[motif == names(motif_library)[i], ])
  plot(log(pval_ref) ~ log_lik_ref, data = motif_pval[motif == names(motif_library)[i], ])
  plot(log(pval_snp) ~ log_lik_snp, data = motif_pval[motif == names(motif_library)[i], ])

  par(mfrow = c(1, 3))
  plot(log(pval_diff) ~ log(pval_rank), data = motif_pval[motif == names(motif_library)[i], ])
  plot(log(pval_rank) ~ log(pval_snp), data = motif_pval[motif == names(motif_library)[i], ])
  plot(log(pval_rank) ~ log(pval_ref), data = motif_pval[motif == names(motif_library)[i], ])
  
  ggplot(aes(x = pval_ref, y = pval_snp, color = pval_rank), data = motif_pval[motif == names(motif_library)[i], ]) + geom_point()

  ggplot(aes(x = pval_ref, y = pval_snp, color = pval_diff), data = motif_pval[motif == names(motif_library)[i], ]) + geom_point()

  system.time(save(motif_library, snpInfo, snp_tbl, motif_scores, file = "~/atsnp_git/atSNP/data/example.rda"))
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  tbl1 <- read.table("~/atsnp_git/data/gwas_snp1.txt", stringsAsFactors = FALSE)
  names(tbl1) <- c("chr", "snp", "snpid")
  tbl1 <- tbl1[grep("rs", tbl1$snpid), ]
  tbl1$snp <- as.integer(tbl1$snp)
  tbl1 <- na.omit(tbl1)
  tbl1$chr <- paste("chr", tbl1$chr, sep = "")
  half.window.size <- 30
  ## ad-hocly remove rows with errors
  
  tbl <- tbl1[!(tbl1$snp - 30) %in% c(82128489,63678456,82727520,63223070,82727627,79864449,79013943,63731181,81244884,81261958,63666419,81714991,78897749,63786954,79632349,78152892,63717525,78655420,78831605,63037654,63778330,63264538,181168108,63696199,63718204,79813518,81287979,81629755,63686837,63712574,78649567,82450909), ]
  tbl$chr[tbl$chr == "chr23"] <- "chrX"
  seqvec <- getSeq(Hsapiens, as.character(tbl$chr), start = tbl$snp - half.window.size, end = tbl$snp + half.window.size, as.character = TRUE)

  codes <- seq(4)
  names(codes) <- c("A", "C", "G", "T")
  sequences <- sapply(seqvec, function(x) codes[strsplit(x, "")[[1]]])
  colnames(sequences) <- tbl$snpid
  rownames(sequences) <- NULL
  sequences <- t(na.omit(t(sequences)))
  transition <- .Call("transition_matrix", sequences, package = "atSNP")
  prior <- apply(transition, 1, sum)
  prior <- prior/sum(prior)
  transition <- transition/apply(transition, 1, sum)
  names(prior) <- colnames(transition) <- rownames(transition) <- c("A", "C", "G", "T")
  save(prior, transition, file = "~/atsnp_git/atSNP/data/default_par.rda")
}

## process the data
data(example)

motif_scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 5)

motif_scores <- MatchSubsequence(motif_scores$snp.tbl, motif_scores$motif.scores, ncores = 3, motif.lib = motif_library)

motif_scores[snpid == "rs2511200" & motif == "ALX3_jolma_DBD_M449", ]

len_seq <- sapply(motif_scores$ref_seq, nchar)
snp_pos <- as.integer(len_seq / 2) + 1

i <- which(motif_scores$snpid == "rs2511200" & motif_scores$motif == "ALX3_jolma_DBD_M449")

test_that("Error: reference bases are not the same as the sequence matrix.", {
  expect_equal(sum(snpInfo$sequence_matrix[31, ] != snpInfo$ref_base), 0)
  expect_equal(sum(snpInfo$sequence_matrix[31, ] == snpInfo$snp_base), 0)
})

test_that("Error: log_lik_ratio is not correct.", {
  expect_equal(motif_scores$log_lik_ref - motif_scores$log_lik_snp, motif_scores$log_lik_ratio)
})

test_that("Error: log likelihoods are not correct.", {

  log_lik <- sapply(seq(nrow(motif_scores)),
                        function(i) {
                          motif_mat <- motif_library[[motif_scores$motif[i]]]
                          bases <- snpInfo$sequence_matrix[motif_scores$ref_start[i]:motif_scores$ref_end[i], motif_scores$snpid[i]]
                          if(motif_scores$ref_strand[i] == "-")
                            bases <- 5 - rev(bases)
                          log(prod(
                                   motif_mat[cbind(seq(nrow(motif_mat)),
                                                   bases)]))
                        })

  expect_equal(log_lik, motif_scores$log_lik_ref)

  snp_mat <- snpInfo$sequence_matrix
  snp_mat[cbind(snp_pos, seq(ncol(snp_mat)))] <- snpInfo$snp_base
  log_lik <- sapply(seq(nrow(motif_scores)),
                    function(i) {
                      motif_mat <- motif_library[[motif_scores$motif[i]]]
                      bases <- snp_mat[motif_scores$snp_start[i]:motif_scores$snp_end[i], motif_scores$snpid[i]]
                      if(motif_scores$snp_strand[i] == "-")
                        bases <- 5 - rev(bases)
                      log(prod(
                               motif_mat[cbind(seq(nrow(motif_mat)),
                                               bases)]))
                    })

  expect_equal(log_lik, motif_scores$log_lik_snp)
})

test_that("Error: log_enhance_odds not correct.", {
  
  len_seq <- sapply(motif_scores$ref_seq, nchar)
  snp_pos <- as.integer(len_seq / 2) + 1

  ## log odds for reduction in binding affinity
  
  pos_in_pwm <- snp_pos - motif_scores$ref_start + 1
  neg_ids <- which(motif_scores$ref_strand == "-")
  pos_in_pwm[neg_ids] <- motif_scores$ref_end[neg_ids]- snp_pos[neg_ids] + 1
  snp_base <- sapply(substr(motif_scores$snp_seq, snp_pos, snp_pos), function(x) which(c("A", "C", "G", "T") == x))
  ref_base <- sapply(substr(motif_scores$ref_seq, snp_pos, snp_pos), function(x) which(c("A", "C", "G", "T") == x))
  snp_base[neg_ids] <- 5 - snp_base[neg_ids]
  ref_base[neg_ids] <- 5 - ref_base[neg_ids]
  my_log_reduce_odds <- sapply(seq(nrow(motif_scores)),
                               function(i)
                               log(motif_library[[motif_scores$motif[i]]][pos_in_pwm[i], ref_base[i]]) -
                               log(motif_library[[motif_scores$motif[i]]][pos_in_pwm[i], snp_base[i]])
                               )

  expect_equal(my_log_reduce_odds, motif_scores$log_reduce_odds)

  ## log odds in enhancing binding affinity
  
  pos_in_pwm <- snp_pos - motif_scores$snp_start + 1
  neg_ids <- which(motif_scores$snp_strand == "-")
  pos_in_pwm[neg_ids] <- motif_scores$snp_end[neg_ids]- snp_pos[neg_ids] + 1
  snp_base <- sapply(substr(motif_scores$snp_seq, snp_pos, snp_pos), function(x) which(c("A", "C", "G", "T") == x))
  ref_base <- sapply(substr(motif_scores$ref_seq, snp_pos, snp_pos), function(x) which(c("A", "C", "G", "T") == x))
  snp_base[neg_ids] <- 5 - snp_base[neg_ids]
  ref_base[neg_ids] <- 5 - ref_base[neg_ids]
  my_log_enhance_odds <- sapply(seq(nrow(motif_scores)),
                                function(i)
                                log(motif_library[[motif_scores$motif[i]]][pos_in_pwm[i], snp_base[i]]) -
                                log(motif_library[[motif_scores$motif[i]]][pos_in_pwm[i], ref_base[i]]) 
                               )

  expect_equal(my_log_enhance_odds, motif_scores$log_enhance_odds)
  

})

test_that("Error: the maximum log likelihood computation is not correct.", {

  snp_mat <- snpInfo$sequence_matrix
  snp_mat[cbind(snp_pos, seq(ncol(snp_mat)))] <- snpInfo$snp_base

  .findMaxLog <- function(seq_vec, pwm) {
    snp_pos <- as.integer(length(seq_vec) / 2) + 1
    start_pos <- snp_pos - nrow(pwm) + 1
    end_pos <- snp_pos
    rev_seq <- 5 - rev(seq_vec)
    
    maxLogProb <- -Inf
    for(i in start_pos : end_pos) {
      LogProb <- log(prod(pwm[cbind(seq(nrow(pwm)),
                                    seq_vec[i - 1 + seq(nrow(pwm))])]))
      if(LogProb > maxLogProb)
        maxLogProb <- LogProb
    }
    for(i in start_pos : end_pos) {
      LogProb <- log(prod(pwm[cbind(seq(nrow(pwm)),
                                    rev_seq[i - 1 + seq(nrow(pwm))])]))
      if(LogProb > maxLogProb)
        maxLogProb <- LogProb
    }
    return(maxLogProb)
  }
  
  ## find the maximum log likelihood on the reference sequence
  my_log_lik_ref <- sapply(seq(nrow(motif_scores)),
                           function(x) {
                             seq_vec<- snpInfo$sequence_matrix[, motif_scores$snpid[x]]
                             pwm <- motif_library[[motif_scores$motif[x]]]
                             return(.findMaxLog(seq_vec, pwm))
                           })

  ## find the maximum log likelihood on the SNP sequence

  my_log_lik_snp <- sapply(seq(nrow(motif_scores)),
                           function(x) {
                             seq_vec<- snp_mat[, motif_scores$snpid[x]]
                             pwm <- motif_library[[motif_scores$motif[x]]]
                             return(.findMaxLog(seq_vec, pwm))
                           })
  
  expect_equal(my_log_lik_ref, motif_scores$log_lik_ref)
  expect_equal(my_log_lik_snp, motif_scores$log_lik_snp)
  
})
