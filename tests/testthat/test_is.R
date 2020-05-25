data(example)

scores <- as.matrix(motif_scores$motif.scores[3:4, 4:5])

test_that("Error: quantile function computing are not equivalent.", {
  for (p in c(0.01, 0.1, 0.5, 0.9, 0.99)) {
    delta <-
      .Call("test_find_percentile", c(scores), p, package = "atSNP")
    delta.r <- -sort(-c(scores))[as.integer(p * length(scores)) + 1]
    expect_equal(delta, delta.r)
  }
})


test_that("Error: utility matrix or starting position distribution is incorrect.",
          {
            artifacts <- gen_test_artifacts()
            for (i in range(10)) {
              motif_len <- nrow(artifacts$pwm)
              theta <- runif(1,
                             min = 0.5,
                             max = 2)
              sample_seq_len <-
                sample(c(1, as.integer(motif_len / 2), motif_len - 1, motif_len * 2), size =
                         1) + motif_len
              delta_r <-
                R_gen_utility_matrix(artifacts$pwm, snpInfo$transition, sample_seq_len, theta)
              delta_cpp <-
                .Call(
                  "test_gen_utility_matrix",
                  artifacts$pwm,
                  snpInfo$transition,
                  sample_seq_len,
                  theta,
                  package = "atSNP"
                )
              prob_start_pos_r <-
                R_gen_prob_start_pos(delta_cpp, motif_len, snpInfo$prior)
              prob_start_pos_cpp <-
                .Call("test_gen_prob_start_pos",
                      delta_cpp,
                      motif_len,
                      snpInfo$prior,
                      package = "atSNP")
              expect_equal(delta_r, delta_cpp, 1e-5)
              expect_equal(prob_start_pos_r, prob_start_pos_cpp, 1e-5)
            }
          })


test_that("Error: the scores for samples are not equivalent.", {
  artifacts <- gen_test_artifacts()
  motif_len <- nrow(artifacts$pwm)
  p <- 0.01
  delta <-
    .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <-
    .Call(
      "test_find_theta",
      artifacts$pwm,
      snpInfo$prior,
      snpInfo$transition,
      delta,
      package = "atSNP"
    )
  ## Use R code to generate a random sample
  for (i in seq(10)) {
    sample <- drawonesample(theta)
    sample_score <-
      .Call(
        "test_compute_sample_score",
        artifacts$pwm,
        sample[seq(2 * motif_len - 1)] - 1,
        sample[motif_len * 2] - 1,
        theta,
        package = "atSNP"
      )
    expect_equal(sample[2 * motif_len + 1], sample_score[4])
    expect_equal(sample[2 * motif_len + 2], sample_score[5])
  }
  ## Use C code to generate a random sample
  for (i in seq(10)) {
    sample_seq_len <-
      sample(c(motif_len - 1, 1, motif_len * 2, as.integer(motif_len / 2)), size =
               1) + motif_len
    delta <- t(artifacts$pwm ^ theta)
    delta <- cbind(matrix(
      sum(snpInfo$prior * delta[, 1]),
      nrow = 4,
      ncol = sample_seq_len - motif_len
    ),
    delta)
    sample_seq <-
      .Call(
        "test_importance_sample",
        delta,
        snpInfo$prior,
        artifacts$snpInfo$transition,
        artifacts$pwm,
        theta,
        package = "atSNP"
      )
    start_pos <- sample_seq[sample_seq_len + 1]
    adj_score <- 0
    for (s in seq(sample_seq_len - motif_len + 1) - 1) {
      adj_score <- adj_score + prod(artifacts$pwm[cbind(seq(motif_len),
                                                        sample_seq[s + seq(motif_len)] + 1)]) ^ theta
    }
    adj_score_cond <-
      prod(artifacts$pwm[cbind(seq(motif_len), sample_seq[start_pos + seq(motif_len)] + 1)]) ^ theta
    sample_score <-
      .Call(
        "test_compute_sample_score",
        artifacts$pwm,
        sample_seq[seq_len(sample_seq_len)],
        sample_seq[sample_seq_len + 1],
        theta,
        package = "atSNP"
      )
    expect_equal(adj_score, sample_score[4])
    expect_equal(adj_score_cond, sample_score[5])
  }
})

test_that("Error: compute the normalizing constant.", {
  ## parameters
  artifacts <- gen_test_artifacts()
  p <- 0.01
  delta <-
    .Call("test_find_percentile", scores, p, package = "atSNP")
  theta <-
    .Call(
      "test_find_theta",
      artifacts$pwm,
      artifacts$snpInfo$prior,
      artifacts$snpInfo$transition,
      delta,
      package = "atSNP"
    )
  ##
  const <-
    .Call(
      "test_func_delta",
      artifacts$pwm,
      artifacts$snpInfo$prior,
      artifacts$snpInfo$transition,
      theta,
      package = "atSNP"
    )
  const.r <-
    sum(
      R_gen_utility_matrix(
        artifacts$pwm,
        artifacts$snpInfo$transition,
        2 * nrow(artifacts$pwm) - 1,
        theta
      )[, seq(nrow(artifacts$pwm))] * artifacts$snpInfo$prior
    )
  expect_equal(abs(const - const.r) / const < 1e-5, TRUE)
})

test_that("Error: the chosen pvalues should have the smaller variance.", {
  artifacts <- gen_test_artifacts()
  motif_len <- nrow(artifacts$pwm)
  .structure <- function(pval_mat) {
    id1 <- apply(pval_mat[, c(2, 4)], 1, which.min)
    return(cbind(pval_mat[, c(1, 3)][cbind(seq_along(id1), id1)],
                 pval_mat[, c(2, 4)][cbind(seq_along(id1), id1)]))
  }
  for (p in c(0.01, 0.05, 0.1)) {
    theta <-
      .Call(
        "test_find_theta",
        artifacts$pwm,
        snpInfo$prior,
        artifacts$snpInfo$transition,
        quantile(c(scores), 1 - p),
        package = "atSNP"
      )
    p_values <-
      .Call(
        "compute_p_values",
        artifacts$pwm,
        snpInfo$prior,
        snpInfo$transition,
        c(scores),
        theta,
        100,
        2 * nrow(artifacts$pwm) - 1,
        0,
        package = "atSNP"
      )
    p_values_s <- .structure(p_values)
    expect_equal(p_values_s[, 2], apply(p_values[, c(2, 4)], 1, min))
  }
})
