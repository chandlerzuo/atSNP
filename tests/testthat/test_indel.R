library(atSNP)
library(BiocParallel)
library(testthat)
data(example)


# Artifacts
trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
test_pwm <- motif_library$SIX5_disc1
adj_pwm <- (test_pwm+0.25) / apply(test_pwm+0.25, 1, sum)
scores <- as.matrix(motif_scores$motif.scores[3:4, 4:5])
score_diff <- abs(scores[, 2] - scores[, 1])


test_importance_sample_indel <- function(
  snpInfo,
  pwm,
  insertion_len,
  score_percentile,
  pwm,
  loglik_type
) {
  adj_pwm <- (pwm+0.25) / apply(pwm+0.25, 1, sum)
  mat_d <- pwm
  for (i in seq(nrow(test_score))) {
    for (j in seq(ncol(test_score))) {
      mat_d[i, j] <- exp(mean(log(pwm[i, j] / pwm[i, -j])))
    }
  }
  return(
    .Call("test_importance_sample_indel",
          snpInfo$prior,
          snpInfo$transition,
          adj_pwm,
          mat_d,
          insertion_len,
          score_percentile,
          pwm,
          loglik_type)
  )
}