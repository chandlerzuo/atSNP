library(MotifAnalysis)
library(testthat)

data(example)
trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$matrix[[1]]
scores <- as.matrix(motif_scores[motif == names(motif_library$matrix)[1], list(log_lik_ref, log_lik_snp)])

## these are functions for this test only
drawonesample <- function(theta) {
  delta <- snpInfo$prior * t(test_pwm ^ theta)
  delta <- delta / rep(apply(delta, 2, sum), each = 4)
  sample <- sample(1:4, 19, replace = TRUE, prob = snpInfo$prior)
  id <- sample(seq(10), 1)
  sample[id : (id + 9)] <- apply(delta, 2, function(x) sample(1:4, 1, prob = x))
  sample <- c(sample,
              id,
              prod(test_pwm[cbind(seq(10),
                                                   sample[id : (id + 9)])]))
  return(sample)
}
jointprob <- function(x) prod(test_pwm[cbind(seq(10), x)])
maxjointprob <- function(x) {
  maxp <- -Inf
  p <- -Inf
  for(i in 1:10) {
    p <- jointprob(x[i:(i+9)])
    if(p > maxp)
      maxp <- p
  }
  for(i in 1:10) {
    p <- jointprob(5 - x[(i+9):i])
    if(p > maxp)
      maxp <- p
  }
  return(maxp)
}
get_freq <- function(sample) {
  ids <- cbind(
               rep(sample[20, ], each = 10) + seq(10),
               rep(seq(1000), each = 10))
  sample_motif <- matrix(sample[ids], nrow = 10) + 1
  emp_freq <- matrix(0, nrow = 10, ncol = 4)
  for(i in seq(10)) {
    for(j in seq(4)) {
      emp_freq[i, j] <- sum(sample_motif[i, ] == j)
    }
  }
  emp_freq <- emp_freq / apply(emp_freq, 1, sum)
  return(emp_freq)
}

test_that("Error: quantile function computing are not equivalent.", {
  for(p in c(1, 10, 50, 90, 99) / 100) {
    delta <- .Call("test_find_percentile", scores, p, package = "MotifAnalysis")
    delta.r <- -sort(-scores)[as.integer(p * length(scores)) + 1]
    expect_equal(delta, delta.r)
  }
})

test_that("Error: the scores for samples are not equivalent.", {
  p <- 0.01
  delta <- .Call("test_find_percentile", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, snpInfo$transition, delta, package = "MotifAnalysis")
  ## Use R code to generate a random sample
  for(i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(19)] - 1, sample[20] - 1, package = "MotifAnalysis")
    expect_equal(sample[21], exp(sample_score[2]))
  }
  ## Use C code to generate a random sample
  for(i in seq(10)) {
    delta <- t(test_pwm ^ theta)
    delta <- cbind(matrix(
                          sum(snpInfo$prior * delta[, 1]),
                          nrow = 4, ncol = 9), delta)
    sample <- .Call("test_importance_sample", delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "MotifAnalysis")
    start_pos <- sample[20]
    adj_score <- prod(test_pwm[cbind(seq(10),
                                     sample[start_pos + seq(10)] + 1)])
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(19)], sample[20], package = "MotifAnalysis")
    expect_equal(adj_score, exp(sample_score[2]))
  }
})

test_that("Error: compute the normalizing constant.", {
  ## parameters
  p <- 0.01
  delta <- .Call("test_find_percentile", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, snpInfo$transition, delta, package = "MotifAnalysis")
  ##
  const <- .Call("test_func_delta", test_pwm, snpInfo$prior, trans_mat, theta, package = "MotifAnalysis")
  const.r <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta, 2, sum))
  expect_equal(abs(const - const.r) / const < 1e-5, TRUE)
})

test_that("Error: sample distributions are not expected.", {
  ## parameters
  p <- 0.1
  delta <- .Call("test_find_percentile", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, trans_mat, delta, package = "MotifAnalysis")
  delta <- t(test_pwm ^ theta)
  delta <- cbind(matrix(
                        sum(snpInfo$prior * delta[, 1]),
                        nrow = 4, ncol = 9), delta)
  registerDoMC(4)
  results <- foreach(i = seq(20)) %dopar% {
    ## generate 1000 samples
    sample <- sapply(seq(1000), function(x)
                     .Call("test_importance_sample",
                           delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "MotifAnalysis"))
    emp_freq1 <- get_freq(sample)
    target_freq <- test_pwm ^ theta * snpInfo$prior
    target_freq <- target_freq / apply(target_freq, 1, sum)
    ## generate samples in R
    sample <- sapply(rep(theta, 1000), drawonesample)
    emp_freq2 <- get_freq(sample[seq(20), ] - 1)
    max(abs(emp_freq1 - target_freq)) > max(abs(emp_freq2 - target_freq))
  }
  sum(unlist(results))
  pbinom(sum(unlist(results)), size = 20, prob = 0.5)
})

## Visual checks
if(FALSE) {

  library(ggplot2)
  
  plot(log(y <- sapply(seq(200) / 100 - 1, function(x)
            .Call("test_func_delta", test_pwm, snpInfo$prior, snpInfo$transition, x, package = "MotifAnalysis"))))

## test the theta

  p_values_1 <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, scores, 0.99, package = "MotifAnalysis")
  p_values_9 <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, scores, 0.1, package = "MotifAnalysis")
  p_values_99 <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, scores, 0.01, package = "MotifAnalysis")
  
  par(mfrow = c(1, 3))
  plot(log(p_values_1[, 1]) ~ scores[, 1])
  plot(log(p_values_9[, 1]) ~ scores[, 1])
  plot(log(p_values_99[, 1]) ~ scores[, 1])
  
  
  p_values_9 <- .Call("test_p_value", test_pwm, snpInfo$prior, trans_mat, scores, 0.1, package = "MotifAnalysis")
  p_values_99 <- .Call("test_p_value", test_pwm, snpInfo$prior, trans_mat, scores, 0.01, package = "MotifAnalysis")
  
  pval_test <- function(x) {
    delta <- .Call("test_find_percentile", scores, x, package = "MotifAnalysis")
    theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, trans_mat, delta, package = "MotifAnalysis")
    const <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta, 2, sum))
    print(const)
    sample <- sapply(rep(theta, 10000), drawonesample)
    pr <- apply(sample[1:19, ], 2, maxjointprob)
    wei <- const / sample[21, ] ^ theta
    print(mean(log(sample[21, ])))
    print(mean(log(pr)))
    print(mean(wei))
    pval <- sapply(scores[,1], function(x) sum(wei[log(pr) > x]) / length(pr))
    return(pval)
  }

  pval_99 <- pval_test(0.01)
  pval_9 <- pval_test(0.1)
  
  rbind(quantile(pval_9, seq(10) / 200),
        quantile(p_values_9[, 1], seq(10) / 200),
        quantile(pval_99, seq(10) / 200),
        quantile(p_values_99[, 1], seq(10) / 200))
  
  plot(log(pval_99), log(p_values_99[, 1]))
  abline(0,1)
  
  plot(log(pval_9), log(p_values_9[, 1]))
  abline(0,1)

}
