#include "ImportanceSampleIndel.h"
#include "MotifScore.h"

/*
Compute the probability that a random sequence can get a score higher than 'score'.
 @arg pwm The position weight matrix, with 4 columns corresponding to A, C, G, T.
 @arg stat_dist A vector of length 4 with stationary distributions of A, C, G, T.
 @arg trans_mat A 4 x 4 transition matrix.
 @arg scores A matrix with 2 columns, with each column corresponding to one allele.
 @arg p The upper percentile of the scores which is used as the mean of the importance
 sampling distribution.
 @arg loglik_type The enum type for max, mean or median.
 @return A matrix with 3 columns. The first two columns are the p-values for the
 log-likelihood scores of each allele. The third column are the p-values for the likelihood ratios.
*/
Rcpp::List p_value_change_indel(
    MarkovChainParam mc_param,
    NumericMatrix mat_d,
    int insertion_len,
    NumericMatrix pwm,
    NumericMatrix adj_pwm,
    NumericVector scores,
    NumericVector pval_ratio,
    double score_percentile,
    int sample_size,
    LoglikType loglik_type)
{
    NumericMatrix p_values(scores.size(), 4);
    NumericVector sample_score(5);

    // find the tilting parameter
    ImportanceSampleIndel sampler(mc_param, adj_pwm, mat_d, insertion_len);
    sampler.initialize(score_percentile);

    double tol = 1e-10;

    rowwise_l1_normalize(pwm, tol);

    for (int i = 0; i < p_values.nrow(); i++)
    {
        for (int j = 0; j < 4; j++)
        {
            p_values(i, j) = 0;
        }
    }
    double mean_score = 0;
    NumericMatrix weights(sample_size, 2), score(sample_size, 2);
    NumericMatrix score_diff(sample_size, 1);
    for (int i = 0; i < sample_size; i++)
    {
        SampleSequence example = sampler.gen_importance_sample();
        ScorePair sample_score_pair = sampler.comp_score_pair(
            pwm,
            example.sequence,
            loglik_type);
        AdjWeights adj_weights = sampler.gen_importance_sample_weights(example.sequence);
        mean_score += sample_score[4];
        // copy the weights and the scores for each allele
        score(i, 0) = sample_score_pair.long_seq_score;
        score(i, 1) = sample_score_pair.short_seq_score;
        score_diff(i, 0) = score(i, 0) - score(i, 1);
        weights(i, 0) = adj_weights.joint;
        weights(i, 1) = adj_weights.base;
    }

    NumericMatrix pval_loglik(scores.size(), 3);
    pval_loglik = comp_empirical_p_values(scores, weights(_, 0), score_diff);

    // compute the sample log ranks
    NumericMatrix pval_ratio_sam(sample_size, 1);
    for (int i = 0; i < sample_size; i++)
    {
        double pval_sam[2];
        for (int j = 0; j < 2; j++)
        {
            pval_sam[j] = 0;
            for (int i1 = 0; i1 < sample_size; i1++)
            {
                if (score(i1, j) >= score(i, j))
                {
                    pval_sam[j] += weights(i1, j);
                }
            }
            if (pval_sam[j] < tol)
            {
                pval_sam[j] = tol;
            }
        }
        pval_ratio_sam(i, 0) = log(pval_sam[0]) - log(pval_sam[1]);
        if (pval_ratio_sam(i, 0) < 0)
        {
            pval_ratio_sam(i, 0) = -pval_ratio_sam(i, 0);
        }
    }

    NumericMatrix pval_rank = comp_empirical_p_values(pval_ratio, weights(_, 0), pval_ratio_sam);

    Rcpp::List ret = Rcpp::List::create(
        Rcpp::Named("score") = pval_loglik,
        Rcpp::Named("rank") = pval_rank);
    return (ret);
}

/* Bisect to find the optimal theta based on a target score.
@arg mat_d The D matrix to induce score change.
@arg insertion_len An integer for the length of insertion.
@arg score The target score.
@return The return value should satisfy: d log(comp_norm_const) / d theta = score
*/
void ImportanceSampleIndel::set_theta(double score)
{
    double tol = 0.01;
    double lower_bound = ImportanceSampleIndel::THETA_MIN;
    double upper_bound = ImportanceSampleIndel::THETA_MAX;
    double mid_point = 0;

    // Outer loop: decrease tol.
    // Inner loop: bisect algorithm.
    while (tol > 1e-4)
    {
        bool lower_sign = this->check_norm_const_diff(lower_bound, tol, score);
        bool upper_sign = this->check_norm_const_diff(upper_bound, tol, score);

        if (!lower_sign && !upper_sign)
        {
            this->theta = upper_bound;
            return;
        }
        else if (lower_sign && upper_sign)
        {
            this->theta = lower_bound;
            return;
        }
        // Bisect. We assume that d log(comp_norm_const) / d theta is increasing with theta.
        if (lower_sign || !upper_sign)
        {
            throw std::range_error("Bisecting error.");
        };
        while (upper_bound - lower_bound > tol)
        {
            mid_point = (lower_bound + upper_bound) * 0.5;
            if (this->check_norm_const_diff(mid_point, tol, score))
            {
                upper_bound = mid_point;
            }
            else
            {
                lower_bound = mid_point;
            }
        }
        tol /= 10;
    }
    this->theta = mid_point;
}

/* Compute the normalization constant.*/
void ImportanceSampleIndel::comp_cond_norm_const()
{
    int motif_len = this->mat_d.nrow();
    NumericMatrix delta(4, motif_len);
    this->cond_norm_const = NumericVector(motif_len + this->insertion_len - 1);

    // Sequence: 0, ..., L-2, [ L-1, ..., L+m-2 ], L+m-1, ..., 2L+m-3
    for (int s = 0; s < motif_len + this->insertion_len - 1; ++s)
    {
        this->cond_norm_const[s] = 0;
        for (int c = motif_len - 1; c <= motif_len + this->insertion_len - 2; ++c)
        {
            if (c - s >= 0 && c - s < motif_len)
            {
                for (int j = 0; j < ImportanceSampleIndel::N_LETTERS; ++j)
                {
                    this->cond_norm_const[s] += exp(log(this->mat_d(c - s, j)) * this->theta);
                }
            }
        }
    }
}

/* Generate a random vector according to the importance sampling distribution.*/
SampleSequence ImportanceSampleIndel::gen_importance_sample()
{
    int motif_len = this->mat_d.nrow();
    int sample_seq_len = 2 * motif_len + this->insertion_len - 2;
    // compute the sampling distribution for each coordinate
    // sample a random vector
    // The longer sequence has length 2*motif_len+insertion_len-2
    // The random vector has one extra element used for the
    // starting position of the matching subsequence.
    RNGScope scope;
    NumericVector rv = runif(sample_seq_len + 1);
    // 1. sample the starting position
    int start_pos = sample_start_position(rv[sample_seq_len]);
    // 2. sample the actual vector
    IntegerVector sample_vec(sample_seq_len);
    for (int i = 0; i < sample_seq_len; i++)
    {
        NumericVector cond_prob(ImportanceSampleIndel::N_LETTERS);
        // 2.1 sample for the stationary distribution part
        if (i < start_pos || i > start_pos + motif_len - 1)
        {
            for (int j = 0; j < ImportanceSampleIndel::N_LETTERS; j++)
            {
                cond_prob[j] = 0;
                if (i == 0 || i == start_pos + motif_len)
                {
                    cond_prob[j] = this->mc_param.stat_dist[j];
                }
                else
                {
                    cond_prob[j] = this->mc_param.trans_mat(sample_vec[i - 1], j);
                }
            }
        }
        // 2.2 Sample for the nonstationary part.
        // Range A: start_pos, ..., start_pos+motif_len-1.
        // Range B: motif_len-1, ..., motif_len+insertion_len-2.
        // For A\B, simulate according to adj_mat.
        // For A intersect B, simulate according to mat_d.
        else
        {
            for (int j = 0; j < ImportanceSampleIndel::N_LETTERS; j++)
            {
                cond_prob[j] = 0;
                if (i >= motif_len - 1 && i <= motif_len + this->insertion_len - 2)
                {
                    if (i - start_pos >= motif_len)
                    {
                        throw std::range_error("Invalid position for the motif length.");
                    };
                    cond_prob[j] = exp(log(this->mat_d(i - start_pos, j)) * this->theta);
                }
                else
                {
                    cond_prob[j] = this->adj_pwm(i - start_pos, j);
                }
            }
        }
        sample_vec[i] = sample_discrete(rv[i], cond_prob);
    }
    SampleSequence seq = {start_pos, sample_vec};
    return (seq);
}

ScorePair ImportanceSampleIndel::comp_score_pair(
    NumericMatrix pwm,
    IntegerVector sample_vec,
    LoglikType loglik_type)
{
    // compute the reverse strand sequence
    int motif_len = pwm.nrow();
    if (pwm.nrow() != ImportanceSampleIndel::N_LETTERS ||
        motif_len != this->mat_d.nrow())
    {
        throw std::length_error("Inconsistent matrix/vector dimensions.");
    };
    // The longer sequence, passed in by 'example', has length
    // motif_len*2+insertion_len-2.
    // The shorter sequence has length motif_len*2-2
    int short_seq_len = motif_len * 2 - 2;
    IntegerVector short_seq(short_seq_len);
    for (int i = 0; i < motif_len - 1; ++i)
    {
        short_seq[i] = sample_vec[i];
        short_seq[short_seq_len - 1 - i] = sample_vec[sample_vec.size() - 1 - i];
    }
    // compute the maximum score for both sequences
    double long_seq_score = 0, short_seq_score = 0;
    SequenceScores long_seq_scores = comp_seq_scores(pwm, sample_vec);
    SequenceScores short_seq_scores = comp_seq_scores(pwm, short_seq);
    switch (loglik_type)
    {
    case LoglikType::mean:
        long_seq_score = long_seq_scores.mean_log_lik;
        short_seq_score = short_seq_scores.mean_log_lik;
        break;
    case LoglikType::median:
        long_seq_score = long_seq_scores.median_log_lik;
        short_seq_score = short_seq_scores.median_log_lik;
        break;
    default:
        long_seq_score = long_seq_scores.max_log_lik;
        short_seq_score = short_seq_scores.max_log_lik;
        break;
    }

    ScorePair sp = {long_seq_score, short_seq_score};
    return (sp);
}

/* Compute the adjustment weights for an importance sample instance.
@arg example The integer vector for the importance sample instance.
*/
AdjWeights ImportanceSampleIndel::gen_importance_sample_weights(
    IntegerVector example)
{
    for (int i = 0; i < example.size(); ++i)
    {
        if (example[i] < 0 || example[i] >= ImportanceSampleIndel::N_LETTERS)
        {
            throw std::range_error("Invalid nucleotide code.");
        };
    }
    if (example.size() != 2 * this->mat_d.nrow() + this->insertion_len - 2)
    {
        throw std::length_error("Invalid sequence length.");
    };
    int motif_len = this->mat_d.nrow();
    double joint_adj = 0;
    for (int s = 0; s < motif_len + this->insertion_len - 1; ++s)
    {
        double adj_s = 0;
        // We need to evaluate the subsequence s, ..., s+motif_len
        // s, ..., s+motif_len-1 is induced by mat_d and adj_pwm
        // s+motif_len is stat_dist for the importance sampling, but
        // should be trans_mat.
        // NOTE: s+motif_len doesn't exist when s=motif_len-1, i.e.
        // the right stationary subsequence does not exist.
        for (int i = s; i <= s + motif_len && i < example.size(); ++i)
        {
            if (i < s + motif_len)
            {
                if (i >= motif_len - 1 && i < motif_len + this->insertion_len - 1)
                {
                    adj_s += this->theta * log(this->mat_d(i - s, example[i]));
                }
                else
                {
                    adj_s += log(this->adj_pwm(i - s, example[i]));
                }
            }
            else
            {
                adj_s += log(this->mc_param.stat_dist[example[i]]);
            }
            if (i == 0)
            {
                adj_s -= log(this->mc_param.stat_dist[example[i]]);
            }
            else
            {
                adj_s -= log(this->mc_param.trans_mat(example[i - 1], example[i]));
            }
        }
        joint_adj += exp(adj_s);
    }
    double base_adj = 0;
    // 0, ..., L-2, [L-1, ..., L+m-2], L+m-1, ..., 2L+m-3
    for (int s = 0; s < motif_len + this->insertion_len - 1; ++s)
    {
        double adj_s = log(this->cond_norm_const[s]);
        for (int i = s; i < s + motif_len; ++i)
        {
            if (i < motif_len - 1 || i >= motif_len + this->insertion_len - 1)
            {
                adj_s += log(this->adj_pwm(i - s, example[i]));
                if (i == 0)
                {
                    adj_s -= log(this->mc_param.stat_dist(example[i]));
                }
                else if (i == motif_len + this->insertion_len - 1)
                {
                    adj_s -= log(this->mc_param.trans_mat(example[motif_len - 2], example[i]));
                }
                else
                {
                    adj_s -= log(this->mc_param.trans_mat(example[i - 1], example[i]));
                }
            }
        }
        // position after the matching subsequence
        if (s + motif_len <= motif_len + this->insertion_len - 1)
        {
            adj_s += log(this->mc_param.stat_dist(example[motif_len + this->insertion_len - 1]));
            adj_s -= log(this->mc_param.trans_mat(example[motif_len - 2], example[motif_len + this->insertion_len - 1]));
        }
        else if (s + motif_len < example.size())
        {
            adj_s += log(this->mc_param.stat_dist(example[s + motif_len]));
            adj_s -= log(this->mc_param.trans_mat(example[s + motif_len - 1], example[s + motif_len]));
        }
        base_adj += exp(adj_s);
    }
    AdjWeights adj_weights = {this->norm_const / joint_adj, this->norm_const / base_adj};
    return adj_weights;
}

/* Helper function to compute the 1st order derivative of log(comp_norm_const)
and compare to the target score
@arg tol Perturbation amount used to compute derivative.
@arg score The target score.
@return true if the derivative >= the target score; false otherwise.
*/
bool ImportanceSampleIndel::check_norm_const_diff(
    double theta,
    double tol,
    double score)
{
    double old_theta = this->theta;
    this->theta = theta + tol / 2;
    this->comp_cond_norm_const();
    this->comp_norm_const();
    double score_1 = log(this->norm_const);
    this->theta = theta - tol / 2;
    this->comp_cond_norm_const();
    this->comp_norm_const();
    double score_0 = log(this->norm_const);
    this->theta = old_theta;
    if (score_1 - score_0 < score * tol)
    {
        return false;
    }
    return true;
}

/* Given a set of parameters,
1. Compute the parameters of the sampler;
2. Generate a random example;
3. Compute the scores and the adjustment weights for this example. 
Return sampling parameters, random example and adjust weights. This 
function is used for unit tests.
*/
SEXP test_importance_sample_indel(
    SEXP _stat_dist,
    SEXP _trans_mat,
    SEXP _adj_pwm,
    SEXP _mat_d,
    SEXP _insertion_len,
    SEXP _score_percentile,
    SEXP _pwm,
    SEXP _loglik_type)
{
    NumericVector stat_dist(_stat_dist);
    NumericMatrix trans_mat(_trans_mat);
    MarkovChainParam mc_param = {stat_dist, trans_mat};
    NumericMatrix adj_pwm(_adj_pwm);
    NumericMatrix mat_d(_mat_d);
    int insertion_len = as<int>(_insertion_len);
    double score_percentile = as<double>(_score_percentile);
    NumericMatrix pwm(_pwm);
    LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));

    ImportanceSampleIndel sampler(mc_param, adj_pwm, mat_d, insertion_len);
    sampler.initialize(score_percentile);
    SampleSequence example = sampler.gen_importance_sample();
    ScorePair score_pair = sampler.comp_score_pair(pwm, example.sequence, loglik_type);
    AdjWeights adj_weights = sampler.gen_importance_sample_weights(example.sequence);
    Rcpp::List ret = Rcpp::List::create(
        Rcpp::Named("start_pos") = example.start_pos,
        Rcpp::Named("sample_sequence") = example.sequence,
        Rcpp::Named("mat_d") = sampler.mat_d,
        Rcpp::Named("theta") = sampler.theta,
        Rcpp::Named("norm_const") = sampler.norm_const,
        Rcpp::Named("cond_norm_const") = sampler.cond_norm_const,
        Rcpp::Named("long_seq_score") = score_pair.long_seq_score,
        Rcpp::Named("short_seq_score") = score_pair.short_seq_score,
        Rcpp::Named("weight_joint") = adj_weights.joint,
        Rcpp::Named("weight_base") = adj_weights.base);
    return (wrap(ret));
}
