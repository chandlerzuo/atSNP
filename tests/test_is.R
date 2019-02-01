library(atSNP)
library(BiocParallel)
library(testthat)
data(example)


trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
id <- 1
test_pwm <- motif_library[[id]]
##test_pwm <- motif_library[["ALX3_jolma_DBD_M449"]]
scores <- as.matrix(motif_scores$motif.scores[motif == names(motif_library)[id], list(log_lik_ref, log_lik_snp)])

motif_len <- nrow(test_pwm)

## these are functions for this test only
drawonesample <- function(theta) {
  delta <- snpInfo$prior * t(test_pwm ^ theta)
  delta <- delta / rep(apply(delta, 2, sum), each = 4)
  sample <- sample(1:4, 2 * motif_len - 1, replace = TRUE, prob = snpInfo$prior)
  id <- sample(seq(motif_len), 1)
  sample[id : (id + motif_len - 1)] <- apply(delta, 2, function(x) sample(1:4, 1, prob = x))
  sc <- s_cond <- 0
  for(s in seq(motif_len)) {
    sc <- sc + prod(test_pwm[cbind(seq(motif_len),
                                  sample[s : (s + motif_len - 1)])]) ^ theta
  }
  s_cond <- prod(test_pwm[cbind(seq(motif_len),
                                sample[id : (id + motif_len - 1)])]) ^ theta
  sample <- c(sample, id, sc, s_cond)
  return(sample)
}
jointprob <- function(x) prod(test_pwm[cbind(seq(motif_len), x)])
maxjointprob <- function(x) {
  maxp <- -Inf
  p <- -Inf
  for(i in 1:motif_len) {
    p <- jointprob(x[i:(i+motif_len - 1)])
    if(p > maxp)
      maxp <- p
  }
  for(i in 1:motif_len) {
    p <- jointprob(5 - x[(i+motif_len - 1):i])
    if(p > maxp)
      maxp <- p
  }
  return(maxp)
}
get_freq <- function(sample) {
  ids <- cbind(
               rep(sample[motif_len * 2, ], each = motif_len) + seq(motif_len),
               rep(seq(100), each = motif_len))
  sample_motif <- matrix(sample[ids], nrow = motif_len) + 1
  emp_freq <- matrix(0, nrow = motif_len, ncol = 4)
  for(i in seq(motif_len)) {
    for(j in seq(4)) {
      emp_freq[i, j] <- sum(sample_motif[i, ] == j)
    }
  }
  emp_freq <- emp_freq / apply(emp_freq, 1, sum)
  return(emp_freq)
}

test_that("Error: quantile function computing are not equivalent.", {
  for(p in c(1, 10, 50, 90, 99) / 100) {
    delta <- .Call("test_find_percentile", c(scores), p, package = "atSNP")
    delta.r <- -sort(-c(scores))[as.integer(p * length(scores)) + 1]
    expect_equal(delta, delta.r)
  }
})

test_that("Error: the scores for samples are not equivalent.", {
  p <- 0.01
  delta <- .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, snpInfo$transition, delta, package = "atSNP")
  ## Use R code to generate a random sample
  for(i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(2 * motif_len - 1)] - 1, sample[motif_len * 2] - 1, theta, package = "atSNP")
    expect_equal(sample[2 * motif_len + 1], sample_score[2])
    expect_equal(sample[2 * motif_len + 2], sample_score[3])
  }
  ## Use C code to generate a random sample
  for(i in seq(10)) {
    delta <- t(test_pwm ^ theta)
    delta <- cbind(matrix(
                          sum(snpInfo$prior * delta[, 1]),
                          nrow = 4, ncol = motif_len - 1), delta)
    sample <- .Call("test_importance_sample", delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "atSNP")
    start_pos <- sample[motif_len * 2]
    adj_score <- 0
    for(s in seq(motif_len) - 1) {
      adj_score <- adj_score + prod(test_pwm[cbind(seq(motif_len),
                                                   sample[s + seq(motif_len)] + 1)]) ^ theta
    }
    adj_score_cond <- prod(test_pwm[cbind(seq(motif_len), sample[start_pos + seq(motif_len)] + 1)]) ^ theta
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(2 * motif_len - 1)], sample[motif_len * 2], theta, package = "atSNP")
    expect_equal(adj_score, sample_score[2])
    expect_equal(adj_score_cond, sample_score[3])
  }
})

test_that("Error: compute the normalizing constant.", {
  ## parameters
  p <- 0.01
  delta <- .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, snpInfo$transition, delta, package = "atSNP")
  ##
  const <- .Call("test_func_delta", test_pwm, snpInfo$prior, trans_mat, theta, package = "atSNP")
  const.r <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta, 2, sum)) * motif_len
  expect_equal(abs(const - const.r) / const < 1e-5, TRUE)
})

test_that("Error: sample distributions are not expected.", {
  ## parameters
  p <- 0.1
  delta <- .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, trans_mat, delta, package = "atSNP")
  delta <- t(test_pwm ^ theta)
  delta <- cbind(matrix(
                        sum(snpInfo$prior * delta[, 1]),
                        nrow = 4, ncol = motif_len - 1), delta)

  results_i <- function(i) {
    ## generate 100 samples
    sample <- sapply(seq(100), function(x)
                     .Call("test_importance_sample",
                           delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "atSNP"))
    emp_freq1 <- get_freq(sample)
    target_freq <- test_pwm ^ theta * snpInfo$prior
    target_freq <- target_freq / apply(target_freq, 1, sum)
    ## generate samples in R
    sample <- sapply(rep(theta, 100), drawonesample)
    emp_freq2 <- get_freq(sample[seq(2 * motif_len), ] - 1)
    max(abs(emp_freq1 - target_freq)) > max(abs(emp_freq2 - target_freq))
  }

  if(Sys.info()[["sysname"]] == "Windows"){
    snow <- SnowParam(workers = 1, type = "SOCK")
    results<-bpmapply(results_i, seq(20), BPPARAM = snow,SIMPLIFY = FALSE)
  }else{
    results<-bpmapply(results_i, seq(20), BPPARAM = MulticoreParam(workers = 1),
                      SIMPLIFY = FALSE)
  }
  
  print(sum(unlist(results)))
  print(pbinom(sum(unlist(results)), size = 20, prob = 0.5))
})

test_that("Error: the chosen pvalues should have the smaller variance.", {
  .structure <- function(pval_mat) {
    id1 <- apply(pval_mat[, c(2, 4)], 1, which.min)
    return(cbind(
                 pval_mat[, c(1, 3)][cbind(seq_along(id1), id1)],
                 pval_mat[, c(2, 4)][cbind(seq_along(id1), id1)])
           )
  }
  for(p in c(0.01, 0.05, 0.1)) {
      theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, trans_mat, quantile(c(scores), 1 - p), package = "atSNP")
      p_values <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, c(scores), theta, 100, package = "atSNP")
    p_values_s <- .structure(p_values)
    expect_equal(p_values_s[, 2], apply(p_values[, c(2, 4)], 1, min))
  }
})
