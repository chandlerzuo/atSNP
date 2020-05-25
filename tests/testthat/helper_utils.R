#' Generate a random sequence following Markov Chain
gen_mc_sequence <- function(prior, transition, sequence_len) {
  dict_size <- length(prior)
  if (dict_size != ncol(transition) ||
      dict_size != nrow(transition)) {
    stop(
      "Incorrect dimension for Markov params, get ",
      length(prior),
      ", (",
      nrow(transition),
      ", ",
      ncol(transition),
      ")."
    )
  }
  sample_seq <- rep(0, sequence_len)
  sample_seq[1] <-
    sample(
      seq_len(dict_size),
      size = 1,
      replace = TRUE,
      prob = prior
    )
  if (sequence_len > 1) {
    for (i in 2:sequence_len) {
      sample_seq[i] <-
        sample(
          seq_len(dict_size),
          size = 1,
          replace = TRUE,
          prob = transition[sample_seq[i - 1], ]
        )
    }
  }
  return(sample_seq)
}


#' Find the log-lik score for the best matching subsequence to an PWM.
R_motif_score_max <- function(sample_seq, pwm) {
  maxlogp <- -Inf
  motif_len <- nrow(pwm)
  for (start_pos in seq(length(sample_seq) - motif_len + 1)) {
    maxlogp <- max(c(
      maxlogp,
      R_motif_score_subseq(sample_seq, pwm, start_pos, FALSE),
      R_motif_score_subseq(sample_seq, pwm, start_pos, TRUE)
    ))
  }
  return(maxlogp)
}


#' Find the log-lik score for a subsequence to an PWM.
R_motif_score_subseq <-
  function(sample_seq, pwm, start_pos, reverse) {
    motif_len <- nrow(pwm)
    if (reverse) {
      ret <-
        sum(log(pwm[cbind(seq(motif_len), (ncol(pwm) + 1 - rev(sample_seq))[
          start_pos:(start_pos + motif_len - 1)])
        ]))
      return(ret)
    }
    ret <-
      sum(log(pwm[cbind(seq(motif_len), sample_seq[start_pos:(start_pos + motif_len - 1)])]))
    return(ret)
  }


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
