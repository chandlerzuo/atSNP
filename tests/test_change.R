library(atSNP)
library(BiocParallel)
library(testthat)
data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
id <- 2
test_pwm <- motif_library[[id]]
scores <- as.matrix(motif_scores$motif.scores[motif == names(motif_library)[id], list(log_lik_ref, log_lik_snp)])
score_diff <- apply(scores, 1, function(x) abs(diff(x)))

pval_a <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, scores, 0.15, 1000)
pval_ratio <- abs(log(pval_a[seq(nrow(scores)),1]) - log(pval_a[seq(nrow(scores)) + nrow(scores), 1]))

test_score <- test_pwm
for(i in seq(nrow(test_score))) {
  for(j in seq(ncol(test_score))) {
    test_score[i, j] <- exp(mean(log(test_pwm[i, j] / test_pwm[i, -j])))
  }
}

adj_mat <- test_pwm + 0.25
motif_len <- nrow(test_pwm)

## these are functions for this test only
drawonesample <- function(theta) {
    prob_start <- rev(apply(test_score ^ theta, 1, sum) / apply(adj_mat, 1, sum))
    id <- sample(seq(motif_len), 1, prob = prob_start)
    sample <- sample(1:4, 2 * motif_len - 1, replace = TRUE, prob = snpInfo$prior)
    delta <- adj_mat
    delta[motif_len - id + 1, ] <- test_score[motif_len - id + 1, ] ^ theta
    sample[id - 1 + seq(motif_len)] <- apply(delta, 1, function(x) sample(seq(4), 1, prob = x))
    ## compute weight
    sc <- 0
    for(s in seq(motif_len)) {
      delta <- adj_mat
      delta[motif_len + 1 - s, ] <- test_score[motif_len + 1 - s, ] ^ theta
      sc <- sc + prod(delta[cbind(seq(motif_len), sample[s - 1 + seq(motif_len)])]) /
        prod(snpInfo$prior[sample[s - 1 + seq(motif_len)]])
    }
    sample <- c(sample, id, sc)
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
  emp_freq <- matrix(0, nrow = 2 * motif_len - 1, ncol = 4)
  for(i in seq(2 * motif_len - 1)) {
    for(j in seq(4)) {
      emp_freq[i, j] <- sum(sample[i, ] == j - 1)
    }
  }
  emp_freq <- emp_freq / apply(emp_freq, 1, sum)
  return(emp_freq)
}

test_that("Error: quantile function computing are not equivalent.", {
  for(p in c(1, 10, 50, 90, 99) / 100) {
    delta <- .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
    delta.r <- as.double(sort(apply(scores, 1, function(x) abs(diff(x))))[as.integer((1 - p) * nrow(scores)) + 1])
    expect_equal(delta, delta.r)
  }
})

test_that("Error: the scores for samples are not equivalent.", {
  p <- 0.1
  delta <- .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
  theta <- .Call("test_find_theta_change", test_score, adj_mat, delta, package = "atSNP")
  ## Use R code to generate a random sample
  for(i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <- .Call("test_compute_sample_score_change", test_pwm, test_score, adj_mat, sample[seq(2 * motif_len - 1)] - 1, snpInfo$prior, trans_mat, sample[2 * motif_len] - 1, theta, package = "atSNP")
    expect_equal(sample[2 * motif_len + 1], sample_score[1])
    sample1 <- sample2 <- sample3 <- sample
    sample1[motif_len] <- seq(4)[-sample[motif_len]][1]
    sample2[motif_len] <- seq(4)[-sample[motif_len]][2]
    sample3[motif_len] <- seq(4)[-sample[motif_len]][3]
    sample_score_r <- log(maxjointprob(sample[seq(2 * motif_len - 1)])) -
      log(c(maxjointprob(sample1[seq(2 * motif_len - 1)]),
            maxjointprob(sample2[seq(2 * motif_len - 1)]),
            maxjointprob(sample3[seq(2 * motif_len - 1)])))
    expect_equal(sample_score_r, sample_score[2:4])
  }
  
  ## Use C code to generate a random sample
  for(i in seq(10)) {
    sample <- .Call("test_importance_sample_change", test_score, snpInfo$prior, trans_mat, test_pwm, theta, package = "atSNP")
    start_pos <- sample[2 * motif_len] + 1
    adj_score <- 0
    for(s in seq_len(motif_len)) {
      adj_s <- sum(log(adj_mat[cbind(seq(motif_len), sample[s - 1 + seq(motif_len)] + 1)]) -
                   log(snpInfo$prior[sample[s - 1 + seq(motif_len)] + 1]))
      adj_s <- adj_s + theta * log(test_score[motif_len + 1 - s, sample[motif_len] + 1]) -
          log(adj_mat[motif_len + 1 - s, sample[motif_len] + 1])
      adj_score <- adj_score + exp(adj_s)
    }
    sample_score <- .Call("test_compute_sample_score_change", test_pwm, test_score, adj_mat, sample[seq(2 * motif_len - 1)], snpInfo$prior, trans_mat, sample[2 * motif_len], theta, package = "atSNP")
    expect_equal(adj_score, sample_score[1])
  }
})

test_that("Error: compute the normalizing constant.", {
  ## parameters
  for(p in seq(9) / 10) {
    delta <- .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
    theta <- .Call("test_find_theta_change", test_score, adj_mat, delta, package = "atSNP")
    const <- .Call("test_func_delta_change", test_score, adj_mat, theta, package = "atSNP")
    ## in R
    adj_sum <- apply(adj_mat, 1, sum)
    wei_sum <- apply(test_score ^ theta, 1, sum)
    const.r <- prod(adj_sum) * sum(wei_sum / adj_sum)
    expect_equal(const, const.r)
  }
})

test_that("Error: sample distributions are not expected.", {
  ## parameters
  p <- 0.1
  delta <- .Call("test_find_percentile_change", score_diff, p, package = "atSNP")
  theta <- .Call("test_find_theta_change", test_score, adj_mat, delta, package = "atSNP")
  prob_start <- rev(apply(test_score ^ theta, 1, sum) / apply(adj_mat, 1, sum))
  ## construct the delta matrix
  delta <- matrix(1, nrow = 4 * motif_len, ncol = 2 * motif_len - 1)
  for(pos in seq(motif_len)) {
    delta[seq(4) + 4 * (pos - 1), ] <- snpInfo$prior
    delta[seq(4) + 4 * (pos - 1), pos - 1 + seq(motif_len)] <- t(test_pwm)
    delta[seq(4) + 4 * (pos - 1), motif_len] <- test_score[motif_len + 1 - pos, ] ^ theta
    delta[seq(4) + 4 * (pos - 1), ] <- delta[seq(4) + 4 * (pos - 1),] / rep(apply(delta[seq(4) + 4 * (pos - 1), ], 2, sum), each = 4)
  }
  target_freq <- matrix(0, nrow = 4, ncol = 2 * motif_len - 1)
  for(pos in seq(motif_len)) {
    target_freq <- target_freq + delta[seq(4) + 4 * (pos - 1), ] * prob_start[pos]
  }
  target_freq <- t(target_freq)
  target_freq <- target_freq / apply(target_freq, 1, sum)

  results_i <- function(i) {
    ## generate 1000 samples
    sample1 <- sapply(seq(1000), function(x)
      .Call("test_importance_sample_change",
            adj_mat, snpInfo$prior, trans_mat, test_score, theta, package = "atSNP"))
    emp_freq1 <- get_freq(sample1)
    sample2 <- sapply(rep(theta, 1000), drawonesample)
    emp_freq2 <- get_freq(sample2 - 1)
    ##    print(rbind(emp_freq1[10, ], emp_freq2[10, ], target_freq[10, ]))
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
  .structure_diff <- function(pval_mat) {
    id <- apply(pval_mat[, c(2, 4)], 1, which.min)
    return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id), id)],
                 pval_mat[, c(2, 4)][cbind(seq_along(id), id)]))
  }
  for(p in c(0.05, 0.1, 0.2, 0.5)) {
    p_values <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, score_diff, pval_ratio, quantile(score_diff, 1 - p), 1000, package = "atSNP")$score
    p_values_s <- .structure_diff(p_values)
    expect_equal(p_values_s[, 2], apply(p_values[, c(2, 4)], 1, min))
  }
})
         
## Visual checks
if(FALSE) {

  plot(diff(log(y <- sapply(seq(100) / 10000 - 0.035, function(x)
            .Call("test_func_delta_change", test_score, adj_mat, x, package = "atSNP"))))* 10000)

  func_delta <- function(theta) {
      log(.Call("test_func_delta_change", test_score, adj_mat, theta, package = "atSNP"))
  }

  1e4 * (func_delta(-0.165 + 1e-4 / 2) - func_delta(-0.165 - 1e-4 / 2))

  ## test the theta
  p_values_9 <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, score_diff, pval_ratio, quantile(score_diff, 0.9), 1000, package = "atSNP")
  p_values_8 <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, score_diff, pval_ratio, quantile(score_diff, 0.8), 1000, package = "atSNP")

  p_values_1 <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, snpInfo$transition, score_diff, pval_ratio, quantile(score_diff, 0.1), 1000, package = "atSNP")

  test <- .Call("test_p_value_change", test_pwm, test_score, (adj_mat +0.25) / 2, snpInfo$prior, snpInfo$transition, score_diff, pval_ratio, 0, 1000, package = "atSNP")

  plot(p_values_9$score[, 1], p_values_9$score[, 2])

  plot(p_values_9$score[, 1], p_values_9$score[, 3])
  
  plot(p_values_9$score[, 2], p_values_9$score[, 4], xlim = range(p_values_9$score[, c(2, 4)]), ylim = range(p_values_9$score[, c(2, 4)]))
  abline(0, 1)
  
  plot(p_values_8$score[, 2], p_values_8$score[, 4], xlim = range(p_values_8$score[, c(2, 4)]), ylim = range(p_values_8$score[, c(2, 4)]))
  abline(0, 1)

  plot(p_values_8$score[, 2], p_values_9$score[, 2])
  abline(0, 1)
  
  p_values <- cbind(p_values_9$score[, 1], p_values_8$score[, 1], p_values_9$score[, 3], p_values_8$score[, 3])[cbind(seq(nrow(p_values_9$score)), apply(cbind(p_values_9$score[, 2], p_values_8$score[, 2], p_values_9$score[, 4], p_values_8$score[, 4]), 1, which.min))]
  
  par(mfrow = c(1, 3))
  plot(log(p_values_9$score[, 1]) ~ score_diff, ylim = c(-5, 0))
  plot(log(p_values_8$score[, 1]) ~ score_diff, ylim = c(-5, 0))
  plot(log(p_values) ~ score_diff, ylim = c(-5, 0))

  p_values <- cbind(p_values_9$rank[, 1], p_values_8$rank[, 1], p_values_9$rank[, 3], p_values_8$rank[, 3])[cbind(seq(nrow(p_values_9$rank)), apply(cbind(p_values_9$rank[, 2], p_values_8$rank[, 2], p_values_9$rank[, 4], p_values_8$rank[, 4]), 1, which.min))]

  par(mfrow = c(1, 3))
  plot(log(p_values_9$rank[, 1]) ~ pval_ratio, ylim = c(-5, 0))
  plot(log(p_values_8$rank[, 1]) ~ pval_ratio, ylim = c(-5, 0))
  plot(log(p_values) ~ pval_ratio, ylim = c(-5, 0))

  par(mfrow = c(1, 2))
  plot(log(p_values_9$rank[, 1]) ~ log(p_values_9$score[, 1]))
  plot(score_diff ~ pval_ratio)
  
  p_values_9 <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, trans_mat, score_diff, quantile(score_diff, 0.9), 1000, package = "atSNP")
  p_values_8 <- .Call("test_p_value_change", test_pwm, test_score, adj_mat, snpInfo$prior, trans_mat, score_diff, quantile(score_diff, 0.8), 1000, package = "atSNP")
  
  pval_test <- function(x) {
      delta <- .Call("test_find_percentile_change", score_diff, x, package = "atSNP")
      theta <- .Call("test_find_theta_change", test_score, adj_mat, delta, package = "atSNP")
      const <- .Call("test_func_delta_change", test_score, adj_mat, theta, package = "atSNP")
      message("Constant value: ", const)
      nrep <- 1000
      log_diff <- rep(0, 3 * nrep)
      wei <- rep(0, nrep)
      sc <- rep(0, nrep)
##      set.seed(0)
      for(i in seq(nrep)) {
          sample <- drawonesample(theta)
          sample_score <- .Call("test_compute_sample_score_change", test_pwm, test_score, adj_mat, sample[seq(2 * motif_len - 1)] - 1, snpInfo$prior, trans_mat, sample[2 * motif_len] - 1, theta, package = "atSNP")

          sample1 <- sample2 <- sample3 <- sample
          sample1[motif_len] <- seq(4)[-sample[motif_len]][1]
          sample2[motif_len] <- seq(4)[-sample[motif_len]][2]
          sample3[motif_len] <- seq(4)[-sample[motif_len]][3]
          pr1 <- maxjointprob(sample1[seq(2 * motif_len - 1)])
          pr2 <- maxjointprob(sample2[seq(2 * motif_len - 1)])
          pr3 <- maxjointprob(sample3[seq(2 * motif_len - 1)])
          pr <- maxjointprob(sample[seq(2 * motif_len - 1)])
          sample_score_r <- c(sample[2 * motif_len + 1], log(pr) - log(c(pr1, pr2, pr3)))
          expect_equal(sample_score[1:4], sample_score_r)
          ## if use sample_score[-1], the result is the same as .Call
          ## if use sample_score_r[-1], the result is the same as pval_test
          log_diff[seq(3) + 3 * (i - 1)] <- sample_score[2:4]
          wei[i] <- const / sample_score[1]
          sc[i] <- log(test_score[motif_len + 1 - sample[2 * motif_len], sample[motif_len]])
      }
      message("Mean weight: ", mean(wei))
      message("Mean diff score: ", mean(log_diff))
      pval <- sapply(score_diff, function(x) sum(rep(wei, each = 3)[abs(log_diff) >= x]) / length(log_diff))
      message("Mean tilting score: ", mean(sc))
      return(pval)
  }

  pval_8 <- pval_test(0.2)
  
  pval_9 <- pval_test(0.1)

  par(mfrow = c(1, 2))
  plot(log(pval_8), log(p_values_8$score[, 1]))
  abline(0,1)  
  plot(log(pval_9), log(p_values_9$score[, 1]))
  abline(0,1)

}

