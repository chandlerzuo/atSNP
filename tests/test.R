library(MotifAnalysis)
library(testthat)

if(FALSE) {
  # construct the test data set
  motif_file <- "/p/keles/ENCODE-CHARGE/volume1/ENCODE-Motifs/encode_motifs_for_fimo.txt"
  system.time(motif_library <- LoadMotifLibrary(motif_file))
  system.time(snpInfo <- LoadSNPData("/p/keles/ENCODE-CHARGE/volume2/SNP/hg19_allinfo.bed", nrow = 100))
  motif_library$matrix <- motif_library$matrix[1:5]
  motif_scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 5)
  system.time(save(motif_library, snpInfo, motif_scores, file = "~/MotifAnalysis_git/MotifAnalysis/data/example.rda"))
}

## process the data
data(example)

motif_scores <- ComputeMotifScore(motif_library, snpInfo, ncores = 5)

len_seq <- sapply(motif_scores$ref_seq, nchar)
snp_pos <- as.integer(len_seq / 2) + 1

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
                          motif_mat <- motif_library$matrix[[motif_scores$motif[i]]]
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
                      motif_mat <- motif_library$matrix[[motif_scores$motif[i]]]
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
                               log(motif_library$matrix[[motif_scores$motif[i]]][pos_in_pwm[i], ref_base[i]]) -
                               log(motif_library$matrix[[motif_scores$motif[i]]][pos_in_pwm[i], snp_base[i]])
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
                                log(motif_library$matrix[[motif_scores$motif[i]]][pos_in_pwm[i], snp_base[i]]) -
                                log(motif_library$matrix[[motif_scores$motif[i]]][pos_in_pwm[i], ref_base[i]]) 
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
                             pwm <- motif_library$matrix[[motif_scores$motif[x]]]
                             return(.findMaxLog(seq_vec, pwm))
                           })

  ## find the maximum log likelihood on the SNP sequence

  my_log_lik_snp <- sapply(seq(nrow(motif_scores)),
                           function(x) {
                             seq_vec<- snp_mat[, motif_scores$snpid[x]]
                             pwm <- motif_library$matrix[[motif_scores$motif[x]]]
                             return(.findMaxLog(seq_vec, pwm))
                           })
  
  expect_equal(my_log_lik_ref, motif_scores$log_lik_ref)
  expect_equal(my_log_lik_snp, motif_scores$log_lik_snp)
  
})
