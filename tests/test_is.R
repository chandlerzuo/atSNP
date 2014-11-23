library(MotifAnalysis)

data(example)

trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
id <- 7
test_pwm <- motif_library$matrix[[id]]
##test_pwm <- motif_library$matrix[["ALX3_jolma_DBD_M449"]]
scores <- as.matrix(motif_scores$motif.scores[motif == names(motif_library$matrix)[id], list(log_lik_ref, log_lik_snp)])

motif_len <- nrow(test_pwm)

## these are functions for this test only
drawonesample <- function(theta) {
  delta <- snpInfo$prior * t(test_pwm ^ theta)
  delta <- delta / rep(apply(delta, 2, sum), each = 4)
  sample <- sample(1:4, 2 * motif_len - 1, replace = TRUE, prob = snpInfo$prior)
  id <- sample(seq(motif_len), 1)
  sample[id : (id + motif_len - 1)] <- apply(delta, 2, function(x) sample(1:4, 1, prob = x))
  sample <- c(sample,
              id,
              prod(test_pwm[cbind(seq(motif_len),
                                                   sample[id : (id + motif_len - 1)])]))
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
               rep(seq(1000), each = motif_len))
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
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(2 * motif_len - 1)] - 1, sample[motif_len * 2] - 1, package = "MotifAnalysis")
    expect_equal(sample[2 * motif_len + 1], exp(sample_score[2]))
  }
  ## Use C code to generate a random sample
  for(i in seq(10)) {
    delta <- t(test_pwm ^ theta)
    delta <- cbind(matrix(
                          sum(snpInfo$prior * delta[, 1]),
                          nrow = 4, ncol = motif_len - 1), delta)
    sample <- .Call("test_importance_sample", delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "MotifAnalysis")
    start_pos <- sample[motif_len * 2]
    adj_score <- prod(test_pwm[cbind(seq(motif_len),
                                     sample[start_pos + seq(motif_len)] + 1)])
    sample_score <- .Call("test_compute_sample_score", test_pwm, sample[seq(2 * motif_len - 1)], sample[motif_len * 2], package = "MotifAnalysis")
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
                        nrow = 4, ncol = motif_len - 1), delta)
  registerDoMC(4)
  results <- foreach(i = seq(motif_len * 2)) %dopar% {
    ## generate 1000 samples
    sample <- sapply(seq(1000), function(x)
                     .Call("test_importance_sample",
                           delta, snpInfo$prior, trans_mat, test_pwm, theta, package = "MotifAnalysis"))
    emp_freq1 <- get_freq(sample)
    target_freq <- test_pwm ^ theta * snpInfo$prior
    target_freq <- target_freq / apply(target_freq, 1, sum)
    ## generate samples in R
    sample <- sapply(rep(theta, 1000), drawonesample)
    emp_freq2 <- get_freq(sample[seq(2 * motif_len), ] - 1)
    max(abs(emp_freq1 - target_freq)) > max(abs(emp_freq2 - target_freq))
  }
  print(sum(unlist(results)))
  print(pbinom(sum(unlist(results)), size = 20, prob = 0.5))
})

test_that("Error: the chosen pvalues should have the smaller variance.", {

  .structure <- function(pval_mat) {
    id1 <- apply(pval_mat[, c(2, 4)], 1, which.min)
    id2 <- apply(pval_mat[, c(6, 8)], 1, which.min)
    return(cbind(
                 pval_mat[, c(1, 3)][cbind(seq_along(id1), id1)],
                 pval_mat[, c(5, 7)][cbind(seq_along(id2), id2)],
                 pval_mat[, c(2, 4)][cbind(seq_along(id1), id1)],
                 pval_mat[, c(6, 8)][cbind(seq_along(id2), id2)])
           )
  }

  for(p in c(0.01, 0.05, 0.1)) {
    p_values <- .Call("test_p_value", test_pwm, snpInfo$prior, snpInfo$transition, scores, p, package = "MotifAnalysis")
    p_values_s <- .structure(p_values)
    expect_equal(p_values_s[, 3], apply(p_values[, c(2, 4)], 1, min))
    expect_equal(p_values_s[, 4], apply(p_values[, c(6, 8)], 1, min))
  }
  
})

## Visual checks
if(FALSE) {

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

  par(mfrow = c(1, 3))
  plot(log(p_values_1[, 3]) ~ scores[, 1])
  plot(log(p_values_9[, 3]) ~ scores[, 1])
  plot(log(p_values_99[, 3]) ~ scores[, 1])

  p_values_9 <- .Call("test_p_value", test_pwm, snpInfo$prior, trans_mat, scores, 0.1, package = "MotifAnalysis")
  p_values_99 <- .Call("test_p_value", test_pwm, snpInfo$prior, trans_mat, scores, 0.01, package = "MotifAnalysis")
  
  pval_test <- function(x) {
    delta <- .Call("test_find_percentile", scores, x, package = "MotifAnalysis")
    theta <- .Call("test_find_theta", test_pwm, snpInfo$prior, trans_mat, delta, package = "MotifAnalysis")
    const <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta, 2, sum))
    print(const)
    sample <- sapply(rep(theta, 10000), drawonesample)
    pr <- apply(sample[seq(2 * motif_len - 1), ], 2, maxjointprob)
    wei <- const / sample[2 * motif_len + 1, ] ^ theta
    print(mean(log(sample[2 * motif_len + 1, ])))
    print(mean(log(pr)))
    print(mean(wei))
    pval <- sapply(scores[,1], function(x) sum(wei[log(pr) >= x]) / length(pr))
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

  plot(p_values_9[,1], p_values_99[,1], xlim = c(0, 0.001), ylim = c(0, 0.001))
  abline(0, 1)

  plot(pval_9, pval_99, xlim = c(0, 0.001), ylim = c(0, 0.001))
  abline(0, 1)

  plot(p_values_99[, 2], p_values_99[, 4])
  abline(0, 1)

  plot(p_values_99[, 1], p_values_99[, 3])
  abline(0, 1)

  plot(p_values_99[, 1], p_values_99[, 2])

}

## this correspond to setting p = 0.1 and p = 0.01

theta1 <- 0.1
theta2 <- 0.22
const1 <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta1, 2, sum))
const2 <- prod(apply(snpInfo$prior * t(test_pwm) ^ theta2, 2, sum))
sample1 <- sapply(rep(theta1, 10000), drawonesample)
sample2 <- sapply(rep(theta2, 1000), drawonesample)
pr1 <- apply(sample1[seq(2 * motif_len - 1), ], 2, maxjointprob)
wei1 <- const1 / sample1[2 * motif_len + 1, ] ^ theta1
pr2 <- apply(sample2[seq(2 * motif_len - 1), ], 2, maxjointprob)
wei2 <- const2 / sample2[2 * motif_len + 1, ] ^ theta2

par(mfrow = c(1,2))
hist(log(pr1))
hist(log(pr2))

mean(log(pr1))
mean(log(pr2))

x <- max(scores[,1])
sum(wei2[log(pr2) > x]) / length(pr2)
sum((wei2[log(pr2) > x])^2) / length(pr2) - (sum(wei2[log(pr2) > x]) / length(pr2))^2

pval1 <- sapply(scores[,1], function(x) sum(wei1[log(pr1) > x]) / length(pr1))
pval2 <- sapply(scores[,1], function(x) sum(wei2[log(pr2) > x]) / length(pr2))
pval21 <- sapply(scores[,1], function(x) sum((wei1 * wei1)[log(pr1) > x]) / length(pr1))
pval22 <- sapply(scores[,1], function(x) sum((wei2 * wei2)[log(pr2) > x]) / length(pr2))

var1 <- pval21 - pval1 * pval1
var2 <- pval22 - pval2 * pval2

plot(pval1, pval2, xlim = c(0, 0.01), ylim = c(0, 0.01))
abline(0, 1)

plot(var1[pval1 < 0.01], var2[pval1 < 0.01])
abline(0,1)

plot(pval2[pval2 < 0.01], sqrt(var2[pval2 < 0.01]))

plot(pval1[pval2 < 0.01], sqrt(var1[pval2 < 0.01]))

## conclusion: for small p, the variance is too large
