library(atSNP)
library(BiocParallel)
library(testthat)
data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$SIX5_disc1
scores <- as.matrix(motif_scores$motif.scores[3:4, 4:5])

motif_len <- nrow(test_pwm)

## these are functions for this test only
drawonesample <- function(theta) {
  delta <- snpInfo$prior * t(test_pwm ^ theta)
  delta <- delta / rep(colSums(delta), each = 4)
  sample <-
    sample(1:4,
           2 * motif_len - 1,
           replace = TRUE,
           prob = snpInfo$prior)
  id <- sample(seq(motif_len), 1)
  sample[id:(id + motif_len - 1)] <-
    apply(delta, 2, function(x)
      sample(1:4, 1, prob = x))
  sc <- s_cond <- 0
  for (s in seq(motif_len)) {
    sc <- sc + prod(test_pwm[cbind(seq(motif_len),
                                   sample[s:(s + motif_len - 1)])]) ^ theta
  }
  s_cond <- prod(test_pwm[cbind(seq(motif_len),
                                sample[id:(id + motif_len - 1)])]) ^ theta
  sample <- c(sample, id, sc, s_cond)
  return(sample)
}

get_freq <- function(sample) {
  ids <- cbind(rep(sample[motif_len * 2,], each = motif_len) + seq(motif_len),
               rep(seq(100), each = motif_len))
  sample_motif <- matrix(sample[ids], nrow = motif_len) + 1
  emp_freq <- matrix(0, nrow = motif_len, ncol = 4)
  for (i in seq(motif_len)) {
    for (j in seq(4)) {
      emp_freq[i, j] <- sum(sample_motif[i,] == j)
    }
  }
  emp_freq <- emp_freq / rowSums(emp_freq)
  return(emp_freq)
}

if (TRUE) {
  ## parameters
  p <- 0.1
  delta <-
    .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <-
    .Call("test_find_theta",
          test_pwm,
          snpInfo$prior,
          trans_mat,
          delta,
          package = "atSNP")
  delta <- t(test_pwm ^ theta)
  delta <- cbind(matrix(
    sum(snpInfo$prior * delta[, 1]),
    nrow = 4,
    ncol = motif_len - 1
  ), delta)
  
  results_i <- function(i) {
    ## generate 100 samples
    sample <- sapply(seq(100), function(x)
      .Call(
        "test_importance_sample",
        delta,
        snpInfo$prior,
        trans_mat,
        test_pwm,
        theta,
        package = "atSNP"
      ))
    emp_freq1 <- get_freq(sample)
    target_freq <- test_pwm ^ theta * snpInfo$prior
    target_freq <- target_freq / rowSums(target_freq)
    ## generate samples in R
    sample <- sapply(rep(theta, 100), drawonesample)
    emp_freq2 <- get_freq(sample[seq(2 * motif_len),] - 1)
    max(abs(emp_freq1 - target_freq)) > max(abs(emp_freq2 - target_freq))
  }
  
  if (Sys.info()[["sysname"]] == "Windows") {
    snow <- SnowParam(workers = 1, type = "SOCK")
    results <-
      bpmapply(results_i,
               seq(20),
               BPPARAM = snow,
               SIMPLIFY = FALSE)
  } else{
    results <-
      bpmapply(results_i,
               seq(20),
               BPPARAM = MulticoreParam(workers = 1),
               SIMPLIFY = FALSE)
  }
  
  print(sum(unlist(results)))
  print(pbinom(sum(unlist(results)), size = 20, prob = 0.5))
}
