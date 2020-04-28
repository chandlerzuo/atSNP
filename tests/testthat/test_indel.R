test_that("Error: mismatch between theta and the expected score difference.",
          {
            artifacts <- gen_test_artifacts()
            motif_len <- nrow(artifacts$pwm)
            for (insertion_len in c(1,
                                    as.integer(motif_len / 2),
                                    motif_len,
                                    motif_len + 1,
                                    motif_len * 2)) {
              score_percentile <- runif(n = 1, min = 0, max = 10)
              result_cpp <- test_importance_sample_indel(
                snpInfo = artifacts$snpInfo,
                pwm = artifacts$pwm,
                adj_pwm = artifacts$adj_pwm,
                mat_d = artifacts$mat_d,
                insertion_len = insertion_len,
                score_percentile = score_percentile,
                loglik_type = 0
              )
              tol <- 1e-4
              derivative <- (log(
                R_comp_norm_const(
                  mat_d = artifacts$mat_d,
                  insertion_len = insertion_len,
                  theta = result_cpp$theta + tol / 2
                )
              ) - log(
                R_comp_norm_const(
                  mat_d = artifacts$mat_d,
                  insertion_len = insertion_len,
                  theta = result_cpp$theta - tol / 2
                )
              )) / tol
              if (result_cpp$theta < 1 && result_cpp$theta > 0) {
                expect_equal(derivative, score_percentile, tolerance = 0.03)
              } else if (result_cpp$theta == 1) {
                expect_gte(score_percentile, derivative)
              } else {
                expect_lte(score_percentile, derivative)
              }
            }
          })


test_that("Error: normalization constant is not correct.", {
  artifacts <- gen_test_artifacts()
  motif_len <- nrow(artifacts$pwm)
  for (insertion_len in c(1,
                          as.integer(motif_len / 2),
                          motif_len,
                          motif_len + 1,
                          motif_len * 2)) {
    score_percentile <- runif(n = 1, min = 0, max = 10)
    result_cpp <- test_importance_sample_indel(
      snpInfo = artifacts$snpInfo,
      pwm = artifacts$pwm,
      adj_pwm = artifacts$adj_pwm,
      mat_d = artifacts$mat_d,
      insertion_len = insertion_len,
      score_percentile = score_percentile,
      loglik_type = 0
    )
    cond_norm_const <-
      R_comp_cond_norm_const(
        mat_d = artifacts$mat_d,
        insertion_len = insertion_len,
        theta = result_cpp$theta
      )
    expect_equal(cond_norm_const, result_cpp$cond_norm_const, tolerance = 0.01)
    expect_equal(result_cpp$norm_const, sum(result_cpp$norm_const))
  }
})


test_that("Error: importance sample weights are not correct.", {
  artifacts <- gen_test_artifacts()
  motif_len <- nrow(artifacts$pwm)
  for (insertion_len in c(1,
                          as.integer(motif_len / 2),
                          motif_len,
                          motif_len + 1,
                          motif_len * 2)) {
    score_percentile <- runif(n = 1, min = 0, max = 10)
    result_cpp <- test_importance_sample_indel(
      snpInfo = artifacts$snpInfo,
      pwm = artifacts$pwm,
      adj_pwm = artifacts$adj_pwm,
      mat_d = artifacts$mat_d,
      insertion_len = insertion_len,
      score_percentile = score_percentile,
      loglik_type = 0
    )
    result_r <- R_comp_importance_sample_weights(
      adj_pwm = artifacts$adj_pwm,
      mat_d = artifacts$mat_d,
      theta = result_cpp$theta,
      prior = artifacts$snpInfo$prior,
      transition = artifacts$snpInfo$transition,
      sample_seq = result_cpp$sample_sequence + 1
    )
    expect_equal(result_r$weight_joint,
                 result_cpp$weight_joint,
                 tolerance = 0.01)
    expect_equal(result_r$weight_base, result_cpp$weight_base, tolerance =
                   0.01)
  }
})

test_that("Error: max likelihood scores on the sequence pair are not correct.",
          {
            artifacts <- gen_test_artifacts()
            motif_len <- nrow(artifacts$pwm)
            for (insertion_len in c(1,
                                    as.integer(motif_len / 2),
                                    motif_len,
                                    motif_len + 1,
                                    motif_len * 2)) {
              score_percentile <- runif(n = 1, min = 0, max = 10)
              result_cpp <- test_importance_sample_indel(
                snpInfo = artifacts$snpInfo,
                pwm = artifacts$pwm,
                adj_pwm = artifacts$adj_pwm,
                mat_d = artifacts$mat_d,
                insertion_len = insertion_len,
                score_percentile = score_percentile,
                loglik_type = 0
              )
              expect_equal(
                R_motif_score_max(
                  sample_seq = result_cpp$sample_sequence + 1,
                  pwm = artifacts$pwm
                ),
                result_cpp$long_seq_score,
                tolerance = 1e-4
              )
              expect_equal(
                R_motif_score_max(
                  sample_seq = result_cpp$sample_sequence[-(motif_len:(motif_len + insertion_len -
                                                                         1))] + 1,
                  pwm = artifacts$pwm
                ),
                result_cpp$short_seq_score,
                tolerance = 1e-4
              )
            }
          })


test_that("Error: distribution of the motif scores is not correct.",
          {
            artifacts <- gen_test_artifacts()
            motif_len <- nrow(artifacts$pwm)
            sample_size <- 100
            for (insertion_len in c(1,
                                    as.integer(motif_len / 2),
                                    motif_len,
                                    motif_len + 1,
                                    motif_len * 2)) {
              score_percentile <- runif(n = 1, min = 0, max = 10)
              score_pairs <-
                sapply(seq_len(sample_size), function(x) {
                  result_cpp <- test_importance_sample_indel(
                    snpInfo = artifacts$snpInfo,
                    pwm = artifacts$pwm,
                    adj_pwm = artifacts$adj_pwm,
                    mat_d = artifacts$mat_d,
                    insertion_len = insertion_len,
                    score_percentile = score_percentile,
                    loglik_type = 0
                  )
                  R_sample_seq <- R_gen_importance_sample(
                    snpInfo = artifacts$snpInfo,
                    adj_pwm = artifacts$adj_pwm,
                    mat_d = artifacts$mat_d,
                    insertion_len = insertion_len,
                    cond_norm_const = result_cpp$cond_norm_const,
                    theta = result_cpp$theta
                  )
                  R_long_seq_score <-
                    R_motif_score_max(sample_seq = R_sample_seq, pwm = artifacts$pwm)
                  motif_len <- nrow(artifacts$pwm)
                  R_short_seq_score <-
                    R_motif_score_max(sample_seq = R_sample_seq[-(motif_len:(motif_len + insertion_len - 1))], pwm = artifacts$pwm)
                  return(
                    c(
                      result_cpp$long_seq_score,
                      result_cpp$short_seq_score,
                      R_long_seq_score,
                      R_short_seq_score
                    )
                  )
                })
              message("KS test for long sequence scores: ",
                      ks.test(score_pairs[1, ], score_pairs[3,])$p.value)
              message("KS test for short sequence scores: ",
                      ks.test(score_pairs[2, ], score_pairs[4,])$p.value)
              expect_gte(ks.test(score_pairs[1, ], score_pairs[3,])$p.value, 0.05)
              expect_gte(ks.test(score_pairs[2, ], score_pairs[4,])$p.value, 0.05)
            }
          })
