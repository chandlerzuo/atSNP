library(atSNP)
library(BiocParallel)
library(testthat)


#' Generate artifacts for unit tests.
gen_test_artifacts <- function() {
  data(example)
  trans_mat <- matrix(rep(snpInfo$prior, each = 4), nrow = 4)
  test_pwm <- motif_library$SIX5_disc1
  adj_pwm <- (test_pwm + 0.25) / apply(test_pwm + 0.25, 1, sum)
  mat_d <- test_pwm
  for (i in seq(nrow(test_pwm))) {
    for (j in seq(ncol(test_pwm))) {
      mat_d[i, j] <- exp(mean(log(test_pwm[i, j] / snpInfo$prior)))
    }
  }
  return(list(
    snpInfo = snpInfo,
    pwm = test_pwm,
    adj_pwm = adj_pwm,
    mat_d = mat_d
  ))
}


#' A wrapper of the cpp implementation
test_importance_sample_indel <- function(snpInfo,
                                         pwm,
                                         adj_pwm,
                                         mat_d,
                                         insertion_len,
                                         score_percentile,
                                         loglik_type) {
  return(
    .Call(
      "test_importance_sample_indel",
      snpInfo$prior,
      snpInfo$transition,
      adj_pwm,
      mat_d,
      insertion_len,
      score_percentile,
      pwm,
      loglik_type,
      package="atSNP"
    )
  )
}


#' Compute conditional normalization constant
R_comp_cond_norm_const <- function(mat_d, insertion_len, theta) {
  # 1, ..., L-1, [L, ..., L+insertion_len-1], L+insert_len, ..., 2L+insertion_len-2
  motif_len <- nrow(mat_d)
  cond_norm_const <- rep(0, motif_len + insertion_len - 1)
  # if binding start at s, the match subsequence is s, ..., s+L-1, which
  # overlaps the insertion part max(s,L), ..., min(L+insertion_len-1, s+L-1)
  # This corresponds to max(s,L)-s+1, ..., min(insertion_len, s)+L-s in D
  for (s in seq_along(cond_norm_const)) {
    start_idx <- max(c(s, motif_len)) - s + 1
    end_idx <- min(c(insertion_len, s)) + motif_len - s
    cond_norm_const[s] <-
      sum(exp(log(mat_d[start_idx:end_idx, ]) * theta))
  }
  return(cond_norm_const)
}


#' Compute normalization constant
R_comp_norm_const <- function(mat_d, insertion_len, theta) {
  cond_norm_const <-
    R_comp_cond_norm_const(mat_d, insertion_len, theta)
  return(sum(cond_norm_const))
}


#' Compute the importance sample adjustment weights
R_comp_importance_sample_weights <-
  function(adj_pwm,
           mat_d,
           theta,
           prior,
           transition,
           sample_seq) {
    motif_len <- nrow(adj_pwm)
    n_letters <- ncol(adj_pwm)
    insertion_len <- length(sample_seq) - 2 * (motif_len - 1)
    if (nrow(mat_d) != motif_len || ncol(mat_d) != n_letters) {
      stop("Shape of mat_d is invalid: ", dim(mat_d))
    }
    if (length(prior) != n_letters ||
        ncol(transition) != n_letters ||
        nrow(transition) != n_letters) {
      stop("Shape of Markov parameters is invalid: ",
           length(prior),
           " ",
           dim(transition))
    }
    cond_norm_const <-
      R_comp_cond_norm_const(mat_d = mat_d,
                             insertion_len = insertion_len,
                             theta = theta)
    norm_const <- sum(cond_norm_const)
    
    # Adjustment factor for the joint distribution / long sequence
    joint_adj <- 0
    for (start_pos in seq_len(insertion_len + motif_len - 1)) {
      # compute the intersection between s, ..., s+L-1 and insertion part
      intersect_start <- max(c(start_pos, motif_len))
      intersect_end <-
        min(c(start_pos, insertion_len)) + motif_len - 1
      adj_s <-
        sum(log(adj_pwm[cbind(seq_len(motif_len), sample_seq[start_pos:(start_pos + motif_len - 1)])]))
      adj_s <-
        adj_s - sum(log(adj_pwm[cbind((intersect_start:intersect_end) - start_pos + 1,
                                      sample_seq[intersect_start:intersect_end])]))
      adj_s <-
        adj_s + sum(log(mat_d[cbind((intersect_start:intersect_end) - start_pos + 1,
                                    sample_seq[intersect_start:intersect_end])])) * theta
      if (start_pos == 1) {
        adj_s <- adj_s - log(prior[sample_seq[1]])
      } else {
        adj_s <-
          adj_s - log(transition[sample_seq[start_pos - 1], sample_seq[start_pos]])
      }
      adj_s <-
        adj_s - sum(log(transition[cbind(sample_seq[start_pos:(start_pos + motif_len - 2)],
                                         sample_seq[(start_pos + 1):(start_pos + motif_len - 1)])]))
      if (start_pos + motif_len <= length(sample_seq)) {
        adj_s <-
          adj_s + log(prior[sample_seq[start_pos + motif_len]]) - log(transition[sample_seq[start_pos + motif_len - 1], sample_seq[start_pos + motif_len]])
      }
      joint_adj <- joint_adj + exp(adj_s)
    }
    
    # Adjustment factor for the short sequence
    base_adj <- 0
    for (start_pos in seq_len(motif_len + insertion_len - 1)) {
      # normalization constant coming from marginalize the insertion part
      adj_s <- log(cond_norm_const[start_pos])
      
      # subsequence intersecting the left part of the insertion
      if (start_pos <= motif_len - 1) {
        adj_s <-
          adj_s + sum(log(adj_pwm[cbind(seq_len(motif_len - start_pos),
                                        sample_seq[start_pos:(motif_len - 1)])]))
        if (start_pos == 1) {
          adj_s <- adj_s - log(prior[sample_seq[1]])
        } else {
          adj_s <-
            adj_s - log(transition[sample_seq[start_pos - 1], sample_seq[start_pos]])
        }
        if (start_pos + 1 <= motif_len - 1) {
          idx <- (start_pos + 1):(motif_len - 1)
          adj_s <-
            adj_s - sum(log(transition[cbind(sample_seq[idx - 1], sample_seq[idx])]))
        }
      }
      # subsequence intersecting the right part to the insertion
      if (start_pos + motif_len - 1 >= insertion_len + motif_len) {
        idx <- (insertion_len + motif_len):(start_pos + motif_len - 1)
        adj_s <-
          adj_s + sum(log(adj_pwm[cbind(idx - start_pos + 1, sample_seq[idx])]))
        adj_s <-
          adj_s - log(transition[sample_seq[motif_len - 1], sample_seq[motif_len +
                                                                         insertion_len]])
        if (start_pos + motif_len - 1 > motif_len + insertion_len) {
          idx <- (motif_len + insertion_len + 1):(start_pos + motif_len - 1)
          adj_s <-
            adj_s - sum(log(transition[cbind(sample_seq[idx - 1], sample_seq[idx])]))
        }
        # the first position after the matching subsequence
        if (start_pos + motif_len <= length(sample_seq)) {
          adj_s <- adj_s + log(prior[sample_seq[start_pos + motif_len]])
          adj_s <-
            adj_s - log(transition[sample_seq[start_pos + motif_len - 1], sample_seq[start_pos +
                                                                                       motif_len]])
        }
      } else {
        # If the matching sequence does not overlap the right of the insertion,
        # the probability after the insertion needs to be adjusted.
        adj_s <-
          adj_s + log(prior[sample_seq[insertion_len + motif_len]])
        adj_s <-
          adj_s - log(transition[sample_seq[motif_len - 1], sample_seq[motif_len +
                                                                         insertion_len]])
      }
      base_adj <- base_adj + exp(adj_s)
    }
    names(joint_adj) <- names(base_adj) <- NULL
    return(list(
      weight_joint = norm_const / joint_adj,
      weight_base = norm_const / base_adj
    ))
  }


#' Generate a random sample following the importance sampling distribution
R_gen_importance_sample <- function(snpInfo,
                                    adj_pwm,
                                    mat_d,
                                    insertion_len,
                                    cond_norm_const,
                                    theta) {
  motif_len <- nrow(adj_pwm)
  n_letters <- length(snpInfo$prior)
  if (motif_len != nrow(mat_d) ||
      length(cond_norm_const) != motif_len + insertion_len - 1) {
    stop("Incorrect input dimensions for motif lengths.")
  }
  if (n_letters != ncol(mat_d) || n_letters != ncol(adj_pwm)) {
    stop("Incorrect input dimensiions for number of letters.")
  }
  sample_seq <- rep(0, 2 * motif_len + insertion_len - 2)
  start_pos <-
    sample(seq_len(length(cond_norm_const)), size = 1, prob = cond_norm_const)
  # start_pos, ..., start_pos+motif_len-1 is the match subsequence
  # Before and after this matching subsequence, simulate two Markov Chain sequences.
  if (start_pos > 1) {
    sample_seq[seq_len(start_pos - 1)] <-
      gen_mc_sequence(
        prior = snpInfo$prior,
        transition = snpInfo$transition,
        sequence_len = start_pos - 1
      )
  }
  if (start_pos + motif_len < length(sample_seq)) {
    sample_seq[(start_pos + motif_len):length(sample_seq)] <-
      gen_mc_sequence(
        prior = snpInfo$prior,
        transition = snpInfo$transition,
        sequence_len = length(sample_seq) - start_pos - motif_len + 1
      )
  }
  # Within the matching subsequence, sample according to adj_pwm outside
  # motif_len, ..., motif_len+insertion_len-1, and according to mat_d^theta
  # inside this insertion part.
  for (i in start_pos:(start_pos + motif_len - 1)) {
    if (i >= motif_len && i <= motif_len + insertion_len - 1) {
      sample_prob <- exp(log(mat_d[i - start_pos + 1, ]) + theta)
    }
    else {
      sample_prob <- adj_pwm[i - start_pos + 1, ]
    }
    sample_seq[i] <-
      sample(seq_len(n_letters), size = 1, prob = sample_prob)
  }
  # validate results
  if (any(sample_seq > n_letters) || any(sample_seq < 1)) {
    stop("Generated sequence ", sample_seq, " is invalid.")
  }
  return(sample_seq)
}


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
