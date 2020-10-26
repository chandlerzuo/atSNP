#include "ImportanceSampleIndel.h"
#include "MotifScore.h"

/*
Compute likelihood ratio scores for SNPs' effect on motif matching.
@arg _motif_library The list object containing a 'matrix' component, which is a list of position weight matrices.
@arg _indel_info A list object. Each element contains
sequence: A vector for the long sequence;
insertion_len: The interger value for the insertion length, which corresponds to the middle part of the sequence.
@return A list of matrices where each row corresponds to an indel, and each column corresponds to a motif.
*/
RcppExport SEXP comp_indel_motif_scores(
    SEXP _motif_library,
    SEXP _indel_info,
    SEXP _loglik_type)
{
    Rcpp::List pwms(_motif_library);
    Rcpp::List indel_info(_indel_info);
    LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));

    int n_motifs = pwms.size();
    int n_indels = indel_info.size();

    IntegerMatrix match_pos_short(n_indels, n_motifs);
    IntegerMatrix match_pos_long(n_indels, n_motifs);
    NumericMatrix log_lik_ratio(n_indels, n_motifs);
    NumericMatrix log_lik_short(n_indels, n_motifs);
    NumericMatrix log_lik_long(n_indels, n_motifs);

    double tol = 1e-10;

    //for each snp
    for (int indel_id = 0; indel_id < n_indels; indel_id++)
    {
        //construct reverse sequence
        SEXP _single_indel_info(indel_info[indel_id]);
        Rcpp::List single_indel_info(_single_indel_info);
        SEXP _inserted_sequence(single_indel_info["inserted_sequence"]);
        SEXP _insertion_len(single_indel_info["insertion_len"]);
        IntegerVector R_inserted_sequence(_inserted_sequence);
        int insertion_len = as<int>(_insertion_len);
        int seq_len = R_inserted_sequence.size();

        if ((seq_len - insertion_len) % 2 != 0)
        {
            throw std::length_error("Long sequence should have equal length on both sizes of insertion.");
        }
        // change the sequence to be indexed by 0
        NumericVector inserted_sequence(seq_len);
        for (int i = 0; i < inserted_sequence.size(); ++i)
        {
            inserted_sequence[i] = R_inserted_sequence[i] - 1;
        }

        // for each motif
        for (int motif_id = 0; motif_id < n_motifs; motif_id++)
        {
            SEXP _pwm(pwms[motif_id]);
            NumericMatrix pwm(_pwm);
            int motif_len = pwm.nrow();
            rowwise_l1_normalize(pwm, tol);

            if (2 * motif_len + insertion_len - 2 > seq_len)
            {
                throw std::length_error("Inserted sequence does not have enough length.");
            }

            IntegerVector long_seq(2 * motif_len + insertion_len - 2);
            IntegerVector short_seq(2 * motif_len - 2);

            int offset = (seq_len - long_seq.size()) / 2;
            for (int i = 0; i < long_seq.size(); ++i)
            {
                long_seq[i] = inserted_sequence[i + offset];
                if (i < motif_len - 1)
                {
                    short_seq[i] = long_seq[i];
                }
                else if (i >= motif_len - 1 + insertion_len)
                {
                    short_seq[i - insertion_len] = inserted_sequence[i + offset];
                }
            }

            SequenceScores long_seq_scores = comp_seq_scores(pwm, long_seq);
            SequenceScores short_seq_scores = comp_seq_scores(pwm, short_seq);

            if (long_seq_scores.best_match_pos > 0)
            {
                match_pos_long(indel_id, motif_id) = long_seq_scores.best_match_pos + offset;
            }
            else
            {
                match_pos_long(indel_id, motif_id) = long_seq_scores.best_match_pos - offset;
            }
            if (short_seq_scores.best_match_pos > 0)
            {
                match_pos_short(indel_id, motif_id) = short_seq_scores.best_match_pos + offset;
            }
            else
            {
                match_pos_short(indel_id, motif_id) = short_seq_scores.best_match_pos - offset;
            }

            switch (loglik_type)
            {
            case LoglikType::median:
                log_lik_short(indel_id, motif_id) = short_seq_scores.median_log_lik;
                log_lik_long(indel_id, motif_id) = long_seq_scores.median_log_lik;
                log_lik_ratio(indel_id, motif_id) = long_seq_scores.median_log_lik - short_seq_scores.median_log_lik;
                break;
            case LoglikType::mean:
                log_lik_short(indel_id, motif_id) = short_seq_scores.mean_log_lik;
                log_lik_long(indel_id, motif_id) = long_seq_scores.mean_log_lik;
                log_lik_ratio(indel_id, motif_id) = long_seq_scores.mean_log_lik - short_seq_scores.mean_log_lik;
                break;
            default:
                log_lik_short(indel_id, motif_id) = short_seq_scores.max_log_lik;
                log_lik_long(indel_id, motif_id) = long_seq_scores.max_log_lik;
                log_lik_ratio(indel_id, motif_id) = long_seq_scores.max_log_lik - short_seq_scores.max_log_lik;
                break;
            }
        }
    }

    Rcpp::List ret = Rcpp::List::create(
        Rcpp::Named("match_pos_short") = match_pos_short,
        Rcpp::Named("match_pos_long") = match_pos_long,
        Rcpp::Named("log_lik_ratio") = log_lik_ratio,
        Rcpp::Named("log_lik_short") = log_lik_short,
        Rcpp::Named("log_lik_long") = log_lik_long);
    return (wrap(ret));
}

/*
Compute the probability that a random sequence can get a score higher than 'score'.
@arg mc_param The Markov chain parameters.
@arg mat_d The D matrix to induce score change.
@arg insertion_len An integer for the length of insertion.
@arg pwm The position weight matrix, with 4 columns corresponding to A, C, G, T.
@arg scores A matrix with 2 columns, with each column corresponding to one sequence.
@arg pval_ratio A vector for the p-value ratios.
@arg score_percentile The score percentile used to generate the importance sampling distribution.
@arg sample_size Importance sampling sample size.
@arg loglik_type The enum type for max, mean or median.
@return A list containing two elements:
score contains the p-values for score differences;
rank contains the rank p-values.
*/
RcppExport SEXP p_value_change_indel(
    SEXP _trans_mat,
    SEXP _stat_dist,
    SEXP _mat_d,
    SEXP _insertion_len,
    SEXP _pwm,
    SEXP _adj_pwm,
    SEXP _scores,
    SEXP _pval_ratio,
    SEXP _score_percentile,
    SEXP _sample_size,
    SEXP _loglik_type)
{
    NumericMatrix trans_mat(_trans_mat);
    NumericVector stat_dist(_stat_dist);
    MarkovChainParam mc_param = {stat_dist, trans_mat};
    NumericMatrix mat_d(_mat_d);
    int insertion_len = as<int>(_insertion_len);
    NumericMatrix pwm(_pwm);
    NumericMatrix adj_pwm(_adj_pwm);
    NumericVector scores(_scores);
    NumericVector pval_ratio(_pval_ratio);
    double score_percentile = as<double>(_score_percentile);
    int sample_size = as<int>(_sample_size);
    LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));
    ;

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
    double total_weight[2] = {0, 0};
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
        total_weight[0] += adj_weights.joint;
        total_weight[1] += adj_weights.base;
    }

    NumericMatrix pval_loglik = comp_empirical_p_values(scores, weights(_, 0), score_diff, TestType::two_sided);

    // compute the sample log ranks
    NumericMatrix pval_ratio_sam(sample_size, 1);
    for (int i = 0; i < sample_size; i++)
    {
        double pval_sam[2] = {0, 0};
        for (int j = 0; j < 2; j++)
        {
            for (int i1 = 0; i1 < sample_size; i1++)
            {
                if (i1 != i && score(i1, j) >= score(i, j))
                {
                    pval_sam[j] += weights(i1, j);
                }
            }
            if (pval_sam[j] < tol)
            {
                pval_sam[j] = tol;
            }
        }
        pval_ratio_sam(i, 0) = log(pval_sam[0]) - log(pval_sam[1]) - log(total_weight[0] - weights(i, 0)) + log(total_weight[1] - weights[i, 1]);
    }

    NumericMatrix pval_rank = comp_empirical_p_values(pval_ratio, weights(_, 0), pval_ratio_sam, TestType::two_sided);

    Rcpp::List ret = Rcpp::List::create(
        Rcpp::Named("score") = pval_loglik,
        Rcpp::Named("rank") = pval_rank);
    return (wrap(ret));
}

/* Bisect to find the optimal theta based on a target score.
@arg mat_d The D matrix to induce score change.
@arg insertion_len An integer for the length of insertion.
@arg score The target score.
@return The return value should satisfy: d log(comp_norm_const) / d theta = score
*/
void ImportanceSampleIndel::set_theta(double score)
{
    double lower_bound = ImportanceSampleIndel::THETA_MIN;
    double upper_bound = ImportanceSampleIndel::THETA_MAX;
    double mid_point = 0;

    // Outer loop: decrease tol.
    // Inner loop: bisect algorithm.
    bool lower_sign = this->comp_expected_score_diff(lower_bound) > score;
    bool upper_sign = this->comp_expected_score_diff(upper_bound) > score;

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
    // Bisect. We assume that the function is increasing with theta.
    if (lower_sign || !upper_sign)
    {
        throw std::range_error("Bisecting error.");
    };
    while (upper_bound - lower_bound > 1e-10)
    {
        mid_point = (lower_bound + upper_bound) * 0.5;
        if (this->comp_expected_score_diff(mid_point) > score)
        {
            upper_bound = mid_point;
        }
        else
        {
            lower_bound = mid_point;
        }
    }
    this->theta = mid_point;
}

/* Compute the normalization constant.*/
NumericVector ImportanceSampleIndel::_comp_cond_norm_const(double theta)
{
    int motif_len = this->mat_d.nrow();
    NumericVector cond_norm_const(motif_len + this->insertion_len - 1);

    // Sequence: 0, ..., L-2, [ L-1, ..., L+m-2 ], L+m-1, ..., 2L+m-3
    for (int s = 0; s < motif_len + this->insertion_len - 1; ++s)
    {
        cond_norm_const[s] = 1;
        for (int c = motif_len - 1; c <= motif_len + this->insertion_len - 2; ++c)
        {
            if (c - s >= 0 && c - s < motif_len)
            {
                double tmp = 0;
                for (int j = 0; j < ImportanceSampleIndel::N_LETTERS; ++j)
                {
                    tmp += exp(log(this->mat_d(c - s, j)) * theta);
                }
                cond_norm_const[s] *= tmp;
            }
        }
    }
    return cond_norm_const;
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
    if (pwm.ncol() != ImportanceSampleIndel::N_LETTERS ||
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

/* Helper function to compute the expected score difference.
*/
double ImportanceSampleIndel::comp_expected_score_diff(double theta)
{
    int motif_len = this->mat_d.nrow();
    NumericVector cond_norm_const = this->_comp_cond_norm_const(theta);
    double norm_const = this->_comp_norm_const(cond_norm_const);

    // Sequence: 0, ..., L-2, [ L-1, ..., L+m-2 ], L+m-1, ..., 2L+m-3
    double expected_score_diff = 0;
    for (int s = 0; s < motif_len + this->insertion_len - 1; ++s)
    {
        double cond_score_diff = 0;
        for (int c = motif_len - 1; c <= motif_len + this->insertion_len - 2; ++c)
        {
            if (c - s >= 0 && c - s < motif_len)
            {
                double nume = 0, denom = 0;
                for (int j = 0; j < ImportanceSampleIndel::N_LETTERS; ++j)
                {
                    nume += exp(log(this->mat_d(c - s, j)) * theta) * log(this->mat_d(c - s, j));
                    denom += exp(log(this->mat_d(c - s, j)) * theta);
                }
                cond_score_diff += nume / denom;
            }
        }
        expected_score_diff += cond_score_diff * cond_norm_const[s];
    }
    expected_score_diff /= norm_const;
    return expected_score_diff;
}

/* Given a set of parameters,
1. Compute the parameters of the sampler;
2. Generate a random example;
3. Compute the scores and the adjustment weights for this example. 
Return sampling parameters, random example and adjust weights. This 
function is used for unit tests.
*/
RcppExport SEXP test_importance_sample_indel(
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
