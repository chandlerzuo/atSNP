library(atSNP)
library(BiocParallel)
library(testthat)
data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$SIX5_disc1
scores <- as.matrix(motif_scores$motif.scores[3:4, 4:5])
score_diff <- abs(scores[, 2] - scores[, 1])

test_score <- test_pwm
for (i in seq(nrow(test_score))) {
  for (j in seq(ncol(test_score))) {
    test_score[i, j] <- exp(mean(log(test_pwm[i, j] / test_pwm[i,-j])))
  }
}

adj_mat <- test_pwm + 0.25
motif_len <- nrow(test_pwm)

## these are functions for this test only
drawonesample <- function(theta) {
  prob_start <- rev(rowSums(test_score ^ theta) / rowSums(adj_mat))
  id <- sample(seq(motif_len), 1, prob = prob_start)
  sample <-
    sample(1:4,
           2 * motif_len - 1,
           replace = TRUE,
           prob = snpInfo$prior)
  delta <- adj_mat
  delta[motif_len - id + 1,] <-
    test_score[motif_len - id + 1,] ^ theta
  sample[id - 1 + seq(motif_len)] <-
    apply(delta, 1, function(x)
      sample(seq(4), 1, prob = x))
  ## compute weight
  sc <- 0
  for (s in seq(motif_len)) {
    delta <- adj_mat
    delta[motif_len + 1 - s,] <-
      test_score[motif_len + 1 - s,] ^ theta
    sc <-
      sc + prod(delta[cbind(seq(motif_len), sample[s - 1 + seq(motif_len)])]) /
      prod(snpInfo$prior[sample[s - 1 + seq(motif_len)]])
  }
  sample <- c(sample, id, sc)
  return(sample)
}

get_freq <- function(sample) {
  emp_freq <- matrix(0, nrow = 2 * motif_len - 1, ncol = 4)
  for (i in seq(2 * motif_len - 1)) {
    for (j in seq(4)) {
      emp_freq[i, j] <- sum(sample[i,] == j - 1)
    }
  }
  emp_freq <- emp_freq / rowSums(emp_freq)
  return(emp_freq)
}

if (FALSE) {
  ## parameters
  p <- 0.1
  delta <-
    .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
  theta <-
    .Call("test_find_theta_change",
          test_score,
          adj_mat,
          delta,
          package = "atSNP")
  prob_start <- rev(rowSums(test_score ^ theta) / rowSums(adj_mat))
  ## construct the delta matrix
  delta <- matrix(1, nrow = 4 * motif_len, ncol = 2 * motif_len - 1)
  for (pos in seq(motif_len)) {
    delta[seq(4) + 4 * (pos - 1),] <- snpInfo$prior
    delta[seq(4) + 4 * (pos - 1), pos - 1 + seq(motif_len)] <-
      t(test_pwm)
    delta[seq(4) + 4 * (pos - 1), motif_len] <-
      test_score[motif_len + 1 - pos,] ^ theta
    delta[seq(4) + 4 * (pos - 1),] <-
      delta[seq(4) + 4 * (pos - 1), ] / rep(colSums(delta[seq(4) + 4 * (pos - 1),]), each = 4)
  }
  target_freq <- matrix(0, nrow = 4, ncol = 2 * motif_len - 1)
  for (pos in seq(motif_len)) {
    target_freq <-
      target_freq + delta[seq(4) + 4 * (pos - 1),] * prob_start[pos]
  }
  target_freq <- t(target_freq)
  target_freq <- target_freq / rowSums(target_freq)
  
  results_i <- function(i) {
    ## generate 100 samples
    sample1 <- sapply(seq(100), function(x)
      .Call(
        "test_importance_sample_change",
        adj_mat,
        snpInfo$prior,
        trans_mat,
        test_score,
        theta,
        package = "atSNP"
      ))
    emp_freq1 <- get_freq(sample1)
    sample2 <- sapply(rep(theta, 100), drawonesample)
    emp_freq2 <- get_freq(sample2 - 1)
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
}
