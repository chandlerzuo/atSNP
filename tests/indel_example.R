data(example)
trans_mat <- snpInfo$transition
mc_prior <- snpInfo$prior

# Simulate Indel Information Data
n_indels <- 10
all_indel_info <- list()
for (i in seq_len(n_indels)) {
  insertion_len <- runif(n = 1, min = 1, max = 10)
  indel_info[[i]] <- list(
    inserted_sequence = sample(
      seq_len(4),
      size = max(sapply(motif_library, nrow)) * 2 - 2 + motif_len,
      replace = TRUE
    ),
    insertion_len = insertion_len
  )
}

# Compute motif scores
# This function can compute motif scores for multiple Indels and multiple pwms
motif_scores <- .Call("comp_indel_motif_scores",
                      motif_library,
                      indel_info,
                      1,
                      package = "atSNP")

# p-values
p_values <- list()
row_id <- 1
for (indel_info in all_indel_info) {
  for (pwm in pwms) {
    mat_d <- t(t(pwm) / mc_prior)
    mat_d <- mat_d / apply(mat_d, 1, sum)
    scores <- t(c(motif_scores$log_lik_long[i], motif_scores$log_lik_short[i]))
    score_percentile <- quantile(scores[,1], 0.9)
    p_values[[i]] <- .Call("p_value_change_indel",
                           trans_mat,
                           mc_prior,
                           mat_d,
                           indel_info$insertion_len,
                           pwm,
                           (pwm + 0.25) / 2,
                           scores,
                           pval_ratio,
                           score_percentile,
                           100, # importance sampling sample size
                           1 # loglik type option, 1 stands for mean
                           )
  }
}