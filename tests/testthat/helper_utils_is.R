R_gen_utility_matrix <- function(pwm, trans_mat, seq_len, theta) {
  n_letters <- ncol(pwm)
  motif_len <- nrow(pwm)
  delta <- matrix(0, nrow = n_letters, ncol = seq_len)
  delta[, seq_len] <- pwm[motif_len, ] ^ theta
  for (i in seq(seq_len - 1, 1)) {
    delta[, i] <- trans_mat %*% delta[, i + 1]
    if (i - seq_len + motif_len > 0) {
      delta[, i] <- delta[, i] * pwm[i - seq_len + motif_len, ] ^ theta
    }
  }
  return(delta)
}

R_gen_prob_start_pos <- function(delta, motif_len, stat_dist) {
  return(apply(delta[,rev(seq(ncol(delta)-motif_len+1))] * stat_dist, 2, sum))
}

## these are functions for this test only
drawonesample <- function(theta) {
  artifacts <- gen_test_artifacts()
  motif_len <- nrow(artifacts$pwm)
  delta <- snpInfo$prior * t(artifacts$pwm ^ theta)
  delta <- delta / rep(colSums(delta), each = 4)
  sample <-
    sample(1:4,
           2 * motif_len - 1,
           replace = TRUE,
           prob = snpInfo$prior)
  id <- sample(seq(motif_len), 1)
  sample[id:(id + motif_len - 1)] <-
    apply(delta, 2, function(x)
      sample(1:4, 1, prob = x))
  sc <- s_cond <- 0
  for (s in seq(motif_len)) {
    sc <- sc + prod(artifacts$pwm[cbind(seq(motif_len),
                                        sample[s:(s + motif_len - 1)])]) ^ theta
  }
  s_cond <- prod(artifacts$pwm[cbind(seq(motif_len),
                                     sample[id:(id + motif_len - 1)])]) ^ theta
  sample <- c(sample, id, sc, s_cond)
  return(sample)
}
