library(atSNP)
library(BiocParallel)
library(testthat)
data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$SIX5_disc1
scores <- as.matrix(motif_scores$motif.scores[3:4, 4:5])
score_diff <- abs(scores[, 2] - scores[, 1])


pval_a <-
  .Call("compute_p_values",
        test_pwm,
        snpInfo$prior,
        snpInfo$transition,
        scores,
        0.15,
        100,
        0)
pval_ratio <-
  abs(log(pval_a[seq(nrow(scores)), 1]) - log(pval_a[seq(nrow(scores)) + nrow(scores), 1]))

test_score <- test_pwm
for (i in seq(nrow(test_score))) {
  for (j in seq(ncol(test_score))) {
    test_score[i, j] <- exp(mean(log(test_pwm[i, j] / test_pwm[i, -j])))
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
  delta[motif_len - id + 1, ] <-
    test_score[motif_len - id + 1, ] ^ theta
  sample[id - 1 + seq(motif_len)] <-
    apply(delta, 1, function(x)
      sample(seq(4), 1, prob = x))
  ## compute weight
  sc <- 0
  for (s in seq(motif_len)) {
    delta <- adj_mat
    delta[motif_len + 1 - s, ] <-
      test_score[motif_len + 1 - s, ] ^ theta
    sc <-
      sc + prod(delta[cbind(seq(motif_len), sample[s - 1 + seq(motif_len)])]) /
      prod(snpInfo$prior[sample[s - 1 + seq(motif_len)]])
  }
  sample <- c(sample, id, sc)
  return(sample)
}
jointprob <- function(x)
  prod(test_pwm[cbind(seq(motif_len), x)])
maxjointprob <- function(x) {
  maxp <- -Inf
  p <- -Inf
  for (i in 1:motif_len) {
    p <- jointprob(x[i:(i + motif_len - 1)])
    if (p > maxp)
      maxp <- p
  }
  for (i in 1:motif_len) {
    p <- jointprob(5 - x[(i + motif_len - 1):i])
    if (p > maxp)
      maxp <- p
  }
  return(maxp)
}
get_freq <- function(sample) {
  emp_freq <- matrix(0, nrow = 2 * motif_len - 1, ncol = 4)
  for (i in seq(2 * motif_len - 1)) {
    for (j in seq(4)) {
      emp_freq[i, j] <- sum(sample[i, ] == j - 1)
    }
  }
  emp_freq <- emp_freq / rowSums(emp_freq)
  return(emp_freq)
}

test_that("Error: quantile function computing are not equivalent.", {
  for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
    delta <-
      .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
    delta.r <-
      as.double(sort(abs(scores[, 2] - scores[, 1]))[ceiling((1 - p) * (nrow(scores)))])
    expect_equal(delta, delta.r)
  }
})

test_that("Error: the scores for samples are not equivalent.", {
  p <- 0.1
  delta <-
    .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
  theta <-
    .Call("test_find_theta_change",
          test_score,
          adj_mat,
          delta,
          package = "atSNP")
  ## Use R code to generate a random sample
  for (i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <-
      .Call(
        "test_compute_sample_score_change",
        test_pwm,
        test_score,
        adj_mat,
        sample[seq(2 * motif_len - 1)] - 1,
        snpInfo$prior,
        trans_mat,
        sample[2 * motif_len] - 1,
        theta,
        package = "atSNP"
      )
    expect_equal(sample[2 * motif_len + 1], sample_score[1])
    sample1 <- sample2 <- sample3 <- sample
    sample1[motif_len] <- seq(4)[-sample[motif_len]][1]
    sample2[motif_len] <- seq(4)[-sample[motif_len]][2]
    sample3[motif_len] <- seq(4)[-sample[motif_len]][3]
    sample_score_r <-
      log(maxjointprob(sample[seq(2 * motif_len - 1)])) -
      log(c(
        maxjointprob(sample1[seq(2 * motif_len - 1)]),
        maxjointprob(sample2[seq(2 * motif_len - 1)]),
        maxjointprob(sample3[seq(2 * motif_len - 1)])
      ))
    expect_equal(sample_score_r, sample_score[2:4])
  }
  
  ## Use C code to generate a random sample
  for (i in seq(10)) {
    sample <-
      .Call(
        "test_importance_sample_change",
        test_score,
        snpInfo$prior,
        trans_mat,
        test_pwm,
        theta,
        package = "atSNP"
      )
    start_pos <- sample[2 * motif_len] + 1
    adj_score <- 0
    for (s in seq_len(motif_len)) {
      adj_s <-
        sum(log(adj_mat[cbind(seq(motif_len), sample[s - 1 + seq(motif_len)] + 1)]) -
              log(snpInfo$prior[sample[s - 1 + seq(motif_len)] + 1]))
      adj_s <-
        adj_s + theta * log(test_score[motif_len + 1 - s, sample[motif_len] + 1]) -
        log(adj_mat[motif_len + 1 - s, sample[motif_len] + 1])
      adj_score <- adj_score + exp(adj_s)
    }
    sample_score <-
      .Call(
        "test_compute_sample_score_change",
        test_pwm,
        test_score,
        adj_mat,
        sample[seq(2 * motif_len - 1)],
        snpInfo$prior,
        trans_mat,
        sample[2 * motif_len],
        theta,
        package = "atSNP"
      )
    expect_equal(adj_score, sample_score[1])
  }
})

test_that("Error: compute the normalizing constant.", {
  ## parameters
  for (p in seq(9) / 10) {
    delta <-
      .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
    theta <-
      .Call("test_find_theta_change",
            test_score,
            adj_mat,
            delta,
            package = "atSNP")
    const <-
      .Call("test_func_delta_change",
            test_score,
            adj_mat,
            theta,
            package = "atSNP")
    ## in R
    adj_sum <- rowSums(adj_mat)
    wei_sum <- rowSums(test_score ^ theta)
    const.r <- prod(adj_sum) * sum(wei_sum / adj_sum)
    expect_equal(const, const.r)
  }
})

test_that("Error: the chosen pvalues should have the smaller variance.", {
  .structure_diff <- function(pval_mat) {
    id <- apply(pval_mat[, c(2, 4)], 1, which.min)
    return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
                 pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
  }
  for (p in c(0.05, 0.1, 0.2, 0.5)) {
    p_values <-
      .Call(
        "compute_p_value_change",
        test_pwm,
        test_score,
        adj_mat,
        snpInfo$prior,
        snpInfo$transition,
        score_diff,
        pval_ratio,
        quantile(score_diff, 1 - p),
        100,
        0,
        package = "atSNP"
      )$score
    p_values_s <- .structure_diff(p_values)
    expect_equal(p_values_s[, 2], apply(p_values[, c(2, 4)], 1, min))
  }
})