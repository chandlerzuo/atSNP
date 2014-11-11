library(MotifAnalysis)
library(testthat)
data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$matrix[[1]]
scores <- as.matrix(motif_scores[motif == names(motif_library$matrix)[1], list(log_lik_ref, log_lik_snp)])

test_score <- test_pwm
for(i in seq(nrow(test_score))) {
  for(j in seq(ncol(test_score))) {
    test_score[i, j] <- exp(mean(log(test_pwm[i, j] / test_pwm[i, -j])))
  }
}

adj_mat <- test_pwm ^ 0.5

## these are functions for this test only
drawonesample <- function(theta) {
    prob_start <- sapply(seq(10),
                         function(j)
                             sum(snpInfo$prior * test_score[11 - j, ] ^ theta *
                                     adj_mat[11 - j, ]) /
                                         sum(snpInfo$prior * adj_mat[11 - j, ])
                         )
    id <- sample(seq(10), 1, prob = prob_start)
    sample <- sample(1:4, 19, replace = TRUE, prob = snpInfo$prior)
    delta <- adj_mat
    delta[11 - id, ] <- delta[11 - id, ] * test_score[11 - id, ] ^ theta
    sample[id - 1 + seq(10)] <- apply(delta, 1, function(x)
        sample(seq(4), 1, prob = x * snpInfo$prior))
    sample <- c(sample,
                id,
                sum(log(delta[cbind(seq(10), sample[id - 1 + seq(10)])]))
                )
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
  emp_freq <- matrix(0, nrow = 19, ncol = 4)
  for(i in seq(19)) {
    for(j in seq(4)) {
      emp_freq[i, j] <- sum(sample[i, ] == j - 1)
    }
  }
  emp_freq <- emp_freq / apply(emp_freq, 1, sum)
  return(emp_freq)
}

test_that("Error: quantile function computing are not equivalent.", {
  for(p in c(1, 10, 50, 90, 99) / 100) {
    delta <- .Call("test_find_percentile_diff", scores, p, package = "MotifAnalysis")
    delta.r <- as.double(sort(apply(scores, 1, function(x) abs(diff(x))))[as.integer((1 - p) * nrow(scores)) + 1])
    expect_equal(delta, delta.r)
  }
})

test_that("Error: the scores for samples are not equivalent.", {

  p <- 0.1
  delta <- .Call("test_find_percentile_diff", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta_diff", test_score, adj_mat, snpInfo$prior, snpInfo$transition, delta, package = "MotifAnalysis")

  ## Use R code to generate a random sample
  for(i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <- .Call("test_compute_sample_score_diff", test_pwm, test_score, adj_mat, sample[seq(19)] - 1, sample[20] - 1, theta, package = "MotifAnalysis")
    expect_equal(sample[21], sample_score[1])

    sample1 <- sample2 <- sample3 <- sample
    sample1[10] <- seq(4)[-sample[10]][1]
    sample2[10] <- seq(4)[-sample[10]][2]
    sample3[10] <- seq(4)[-sample[10]][3]
    sample_score_r <- log(maxjointprob(sample[seq(19)])) -
      log(c(maxjointprob(sample1[seq(19)]),
            maxjointprob(sample2[seq(19)]),
            maxjointprob(sample3[seq(19)])))
    expect_equal(sample_score_r, sample_score[-1])
  }
  
  ## Use C code to generate a random sample
  delta <- matrix(1, nrow = 40, ncol = 19)
  for(pos in seq(10)) {
      for(j in (pos + 9) : 1) {
          if(j < pos + 9) {
              delta[4 * (pos - 1) + seq(4), j] <- sum(snpInfo$prior * delta[4 * (pos - 1) + seq(4), j + 1])
          }
          if(j >= pos) {
              delta[4 * (pos - 1) + seq(4), j] <- delta[4 * (pos - 1) + seq(4), j] * adj_mat[j - pos + 1, ]
          }
          if(j == 10) {
              delta[4 * (pos - 1) + seq(4), j] <- delta[4 * (pos - 1) + seq(4), j] * test_score[j - pos + 1, ] ^ theta
          }
      }
  }
  for(i in seq(10)) {
    sample <- .Call("test_importance_sample_diff", delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "MotifAnalysis")
    start_pos <- sample[20] + 1
    adj_score <- sum(log(adj_mat[cbind(seq(10), sample[start_pos - 1 + seq(10)] + 1)]))
    adj_score <- adj_score + theta * log(test_score[11 - start_pos, sample[10] + 1])
    sample_score <- .Call("test_compute_sample_score_diff", test_pwm, test_score, adj_mat, sample[seq(19)], sample[20], theta, package = "MotifAnalysis")
    expect_equal(adj_score, sample_score[1])
  }
  
})

test_that("Error: compute the normalizing constant.", {
  
  ## parameters
  p <- 0.1
  delta <- .Call("test_find_percentile_diff", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta_diff", test_score, adj_mat, snpInfo$prior, snpInfo$transition, delta, package = "MotifAnalysis")
  
  ##
  const <- .Call("test_func_delta_diff", test_score, adj_mat, snpInfo$prior, trans_mat, theta, package = "MotifAnalysis")

   prob_start <- sapply(seq(10),
                         function(j)
                             sum(snpInfo$prior * test_score[11 - j, ] ^ theta *
                                     adj_mat[11 - j, ]) /
                                         sum(snpInfo$prior * adj_mat[11 - j, ])
                         )
  
  const.r <- prod(apply(snpInfo$prior * t(adj_mat), 2, sum)) * sum(prob_start)
  expect_equal(const, const.r)

  plot(sapply(seq(-1, 1, by = 0.01), function(x) log(.Call("test_func_delta_diff", test_score, adj_mat, snpInfo$prior, trans_mat, x, package = "MotifAnalysis"))))
})

test_that("Error: sample distributions are not expected.", {
  
  ## parameters
  p <- 0.1
  delta <- .Call("test_find_percentile_diff", scores, p, package = "MotifAnalysis")
  theta <- .Call("test_find_theta_diff", test_score, adj_mat, snpInfo$prior, snpInfo$transition, delta, package = "MotifAnalysis")

  ## construct the delta matrix
  delta <- matrix(1, nrow = 40, ncol = 19)
   for(pos in seq(10)) {
        for(j in (pos + 9) : 1) {
            if(j < pos + 9) {
                delta[4 * (pos - 1) + seq(4), j] <- sum(snpInfo$prior * delta[4 * (pos - 1) + seq(4), j + 1])
            }
            if(j >= pos) {
                delta[4 * (pos - 1) + seq(4), j] <- delta[4 * (pos - 1) + seq(4), j] * adj_mat[j - pos + 1, ]
            }
            if(j == 10) {
                delta[4 * (pos - 1) + seq(4), j] <- delta[4 * (pos - 1) + seq(4), j] * test_score[j - pos + 1, ] ^ theta
            }
        }
    }

  target_freq <- matrix(0, nrow = 4, ncol = 19)
  
  mat <- snpInfo$prior * matrix(delta[, 1], nrow = 4)
  wei <- apply(mat, 2, sum)
  for(j in seq(19)) {
      for(pos in seq(10)) {
          tmp <- delta[seq(4) + 4 * (pos - 1), j] * snpInfo$prior
          target_freq[, j] <- target_freq[, j] +  tmp / sum(tmp) * wei[pos]
      }
  }
  target_freq <- t(target_freq)
  target_freq <- target_freq / apply(target_freq, 1, sum)
  
  registerDoMC(4)
  
  results <- foreach(i = seq(20)) %dopar% {

    ## generate 1000 samples
    sample1 <- sapply(seq(1000), function(x)
                     .Call("test_importance_sample_diff",
                           delta, snpInfo$prior, trans_mat, test_score, theta, package = "MotifAnalysis"))
    emp_freq1 <- get_freq(sample1)
    
    sample2 <- sapply(rep(theta, 1000), drawonesample)
    emp_freq2 <- get_freq(sample2 - 1)

##    print(rbind(emp_freq1[10, ], emp_freq2[10, ], target_freq[10, ]))
    max(abs(emp_freq1 - target_freq)) > max(abs(emp_freq2 - target_freq))
    
  }
  print(sum(unlist(results)))

  print(pbinom(sum(unlist(results)), size = 20, prob = 0.5))
  
})

## Visual checks
if(FALSE) {
  
  plot(log(y <- sapply(seq(200) / 100 - 1, function(x)
            .Call("test_func_delta_diff", test_score, adj_mat, snpInfo$prior, snpInfo$transition, x, package = "MotifAnalysis"))))

## test the theta

  p_values_9 <- .Call("test_p_value_diff", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, scores, 0.1, package = "MotifAnalysis")
  p_values_8 <- .Call("test_p_value_diff", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, scores, 0.2, package = "MotifAnalysis")

  score_diff <- apply(scores, 1, function(x) abs(diff(x)))
  
  par(mfrow = c(1, 2))
  plot(log(p_values_9) ~ score_diff)
  plot(log(p_values_8) ~ score_diff)
  
  p_values_9 <- .Call("test_p_value_diff", test_pwm, test_score, adj_mat, snpInfo$prior, trans_mat, scores, 0.1, package = "MotifAnalysis")
  p_values_8 <- .Call("test_p_value_diff", test_pwm, test_score, adj_mat, snpInfo$prior, trans_mat, scores, 0.2, package = "MotifAnalysis")
  
  pval_test <- function(x) {
      delta <- .Call("test_find_percentile_diff", scores, x, package = "MotifAnalysis")
      theta <- .Call("test_find_theta_diff", test_score, adj_mat, snpInfo$prior, trans_mat, delta, package = "MotifAnalysis")
      const <- .Call("test_func_delta_diff", test_score, adj_mat, snpInfo$prior, trans_mat, theta, package = "MotifAnalysis")
      message("Constant value: ", const)
      log_diff <- rep(0, 3000)
      wei <- rep(0, 1000)
##      set.seed(0)
      for(i in seq(1000)) {
          sample <- drawonesample(theta)
          sample_score <- .Call("test_compute_sample_score_diff", test_pwm, test_score, adj_mat, sample[seq(19)] - 1, sample[20] - 1, theta, package = "MotifAnalysis")

          sample1 <- sample2 <- sample3 <- sample
          sample1[10] <- seq(4)[-sample[10]][1]
          sample2[10] <- seq(4)[-sample[10]][2]
          sample3[10] <- seq(4)[-sample[10]][3]
          pr1 <- maxjointprob(sample1[1:19])
          pr2 <- maxjointprob(sample2[1:19])
          pr3 <- maxjointprob(sample3[1:19])
          pr <- maxjointprob(sample[1:19])
          sample_score_r <- c(sample[21], log(pr) - log(c(pr1, pr2, pr3)))
          expect_equal(sample_score, sample_score_r)
          ## if use sample_score[-1], the result is the same as .Call
          ## if use sample_score_r[-1], the result is the same as pval_test
          log_diff[seq(3) + 3 * (i - 1)] <- sample_score[-1]
          wei[i] <- const / exp(sample_score[1])
      }
      message("Mean weight: ", mean(wei))
      message("Mean diff score: ", mean(log_diff))
      pval <- sapply(score_diff, function(x) sum(rep(wei, each = 3)[log_diff >= x]) / length(log_diff))
      return(pval)
  }


  pval_8 <- pval_test(0.2)
  
  pval_9 <- pval_test(0.1)
  
  plot(log(pval_8), log(p_values_8))
  abline(0,1)
  
  plot(log(pval_9), log(p_values_9))
  abline(0,1)

}
