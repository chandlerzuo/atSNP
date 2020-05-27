library(atSNP)
data(example)
trans_mat <- snpInfo$transition
mc_prior <- snpInfo$prior
LOGLIK_TYPES <- c("max", "mean", "median")
loglik_type <- which(LOGLIK_TYPES == "mean") - 1

# Simulate Indel Information Data
n_indels <- 10
all_indel_info <- list()
for (i in seq_len(n_indels)) {
  insertion_len <- sample(seq(10), size = 1)
  all_indel_info[[i]] <- list(
    inserted_sequence = sample(
      seq_len(4),
      size = max(sapply(motif_library, nrow)) * 2 - 2 + insertion_len,
      replace = TRUE
    ),
    insertion_len = insertion_len
  )
}

# Step 1. Compute motif scores
# This function can compute motif scores for multiple Indels and multiple pwms
# The return value is a list of matrices
# In each matrix, each row corresponds to an Indel, each column corresponds
# to a motif.
motif_scores <- .Call("comp_indel_motif_scores",
                      motif_library,
                      all_indel_info,
                      # select the log-lik type here
                      loglik_type,
                      package = "atSNP")

# Step 2 and 3.
results <- list()
result_id <- 1
for (indel_id in seq_along(all_indel_info)) {
  for (motif_id in seq_along(motif_library)) {

    # Step 2. Compute single allele p-values
    indel_info <- all_indel_info[[indel_id]]
    pwm <- motif_library[[motif_id]]
    # Prerequisite
    # scores should be a matrix of 2 columns, each representing the score
    # for the long and the short sequence.
    # Eventually, scores should be a matrix of multiple rows as well.
    # Each row represents an Indel.
    # These indels should have similar motif scores so that their p-values
    # can be calculated together with shared Monte Carlo samples.
    scores <-
      t(c(motif_scores$log_lik_long[indel_id, motif_id], motif_scores$log_lik_short[indel_id, motif_id]))
    p_value_affinity <- rep(0, 2)
    for (j in seq(2)) {
      if (j == 1) {
        # for long sequence
        sample_seq_len <- 2 * nrow(pwm) - 2 +indel_info$insertion_len
        reference_score <-
          motif_scores$log_lik_long[indel_id, motif_id]
      } else  if (j == 2) {
        # j=2 for short sequence
        sample_seq_len <- 2 * nrow(pwm) - 2
        reference_score <-
          motif_scores$log_lik_long[indel_id, motif_id]
      }
      # Compute theta parameter in importance sampling distribution
      theta <- .Call(
        "test_find_theta",
        pwm,
        snpInfo$prior,
        snpInfo$transition,
        # Importance sample scores will have
        # average value of reference score.
        reference_score,
        sample_seq_len,
        package = "atSNP"
      )
      p_value_affinity[j] <- pval_with_less_var(
        # Importance sampling based p-values can be computed in two forms
        # A) \sum_i weight_i*X_i / N
        # B) \sum_i weight_i*X_i / \sum_i weight_i
        # p_val_with_less_var takes a matrix of 4 columns
        # Column 1 uses form A
        # Column 2 is the estimated variance of Column 1
        # Column 3 uses form B
        # Column 4 is the estimated variance of Column 3.
        # For each row, p_val_with_less_var picks either Column 1 or 3
        # with the smaller estimated variance.
        # The return of p_val_with_less_var is a matrix of 2 column.
        # Column 1 is the estimate in form A or B. Column 2 is the
        # estimated variance.
        .Call(
          "compute_p_values",
          # PWM
          pwm,
          # MC stationary distribution
          snpInfo$prior,
          # transition matrix
          snpInfo$transition,
          # motif score
          # This is an array, representing scores for multiple Indels.
          # In this example, we compute for one Indel at a time.
          # Eventually we should calculate for multiple Indels together,
          # which (1) have similar scores as reference_score;
          # (2) have the same insertion length.
          scores[, j],
          # theta parameter in importance sampling
          theta,
          # Monte-carlo Sample Size
          # We should adaptively use larger sample sizes for large scores.
          # Because those Indels have small p-values, they are more important
          # biologically, and it is better to compute p-values more accurately for them.
          # For small score difference, p-values will be insignificant,
          # so it is better to use small sample size for computation speed.
          # In atSNP, we choose between 2000 and 1e5 based on score differences.
          # Please refer to the implementation for atSNP for details.
          200,
          # The sequence length
          sample_seq_len,
          # Use 1 for mean log lik scores
          loglik_type,
          package = "atSNP"
        )[, seq(4)] # the first 4 columns are p-values
        # the last 4 columns are conditional p-values and are not useful here
      )[, 1]
    }
    
    # Step 3. Compute p-value for motif score change and p-value change.
    # mat_d is the matrix to induce binding affinity change
    mat_d <- t(t(pwm) / mc_prior)
    score_diff <- c(scores[, 1] - scores[, 2])
    # reference_score is used to compute the theta parameter in importance sampling
    reference_score <-
      motif_scores$log_lik_long[indel_id, motif_id] - motif_scores$log_lik_short[indel_id, motif_id]
    p_value_change <-
      .Call(
        "p_value_change_indel",
        # Markov Chain transition matrix
        trans_mat,
        # Markov Chain stationary distribution
        mc_prior,
        # The D matrix used to induce binding affinity change
        mat_d,
        # Insertion length
        indel_info$insertion_len,
        # PWM
        pwm,
        # Adjusted PWM
        (pwm + 0.25) / 2,
        # An array of motif scores.
        # Each entry corresponds to one Indel.
        # Eventually, we should group Indels with similar scores
        # and compute their p-values together.
        score_diff,
        # Difference between p_values
        # This should be an array of the same length as score_diff.
        # Each entry of this array corresponds to one Indel.
        c(p_value_affinity[1] - p_value_affinity[2]),
        # This is used to compute the theta parameter in importance sampling.
        reference_score,
        # importance sampling sample size
        # We should adaptively use larger sample sizes for large score difference,
        # Because those Indels have small p-values, they are more important
        # biologically, and it is better to compute p-values more accurately for them.
        # For small score difference, p-values will be insignificant,
        # so it doesn't matter much if it is not very accurate.
        # In atSNP, we choose between 2000 and 1e5 based on score differences.
        # Please refer to the implementation for atSNP for details.
        100,
        # loglik type option
        loglik_type
      )
    
    # The output format is a list.
    # This can be changed to whatever is convenient for downstream analysis.
    results[[result_id]] <- list(
      motif_scores = scores,
      p_value_change = list(
        rank = pval_with_less_var(p_value_change$rank)[, 1],
        score = pval_with_less_var(p_value_change$score)[, 1]
      ),
      p_value_affinity = p_value_affinity
    )
    result_id <- result_id + 1
  }
}
