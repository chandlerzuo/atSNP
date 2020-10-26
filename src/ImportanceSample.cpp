#include "ImportanceSample.h"

/*
Reference: Hock Peng Chan et al.(2010). Importance sampling of word patterns in DNA and protein sequences. Journal of computational biology, 17(12).
*/

/*
Compute the probability that a random sequence can get a score higher than 'score'.
@arg pwm The position weight matrix, with 4 columns corresponding to A, C, G, T.
@arg stat_dist A vector of length 4 with stationary distributions of A, C, G, T.
@arg trans_mat A 4 x 4 transition matrix.
@arg scores A vector of log-lik scores from both alleles.
@arg theta A float parameter used to construct the importance sampling distribution.
@arg n_sample An integer for the number of Monte Carlo examples.
@arg seq_len An integer for the sequence length. This should be
	SNP: motif length *2 - 1 for SNP;
	Shorter sequence in Indel: 2 * motif length - 2;
	Longer sequence in Indel: 2 * motif length - 2 + insertion_length.
for the longer sequence in Indel.
@arg loglik_type The enum type for max, mean or median.
@return A matrix with 8 columns:
  Column 1-2 are simple estimates of p-values and their variances.
  Column 3-4 are ratio estimates of p-values and their variances.
  Column 5-6 are simple estimates of conditional p-values and their variances.
  Column 7-8 are ratio estimates for conditional p-values and their variances.
*/
NumericMatrix p_value(
	NumericMatrix pwm,
	NumericVector stat_dist,
	NumericMatrix trans_mat,
	NumericVector scores,
	double theta,
	int n_sample,
	int seq_len,
	LoglikType loglik_type)
{
	NumericMatrix p_values(scores.size(), 8);

	double tol = 1e-10;
	int motif_len = pwm.nrow();
	IntegerVector sample_vec(seq_len);
	IntegerVector sample(seq_len + 1);
	NumericVector sample_score(5);

	if (seq_len < motif_len)
	{
		throw std::length_error("Sequence length cannot be smaller than motif length.");
	}

	for (int i = 0; i < 4; i++)
		for (int m = 0; m < motif_len; m++)
			if (pwm(m, i) < tol)
				pwm(m, i) = tol;

	NumericMatrix delta = gen_utility_matrix(pwm, trans_mat, seq_len, theta);
	double norm_const = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < seq_len - motif_len + 1; j++)
		{
			norm_const += stat_dist[i] * delta(i, j);
		}
	}

	for (int i = 0; i < p_values.nrow(); i++)
		for (int j = 0; j < 4; j++)
			p_values(i, j) = 0;
	double mean_sample = 0;
	double mean_adj_score = 0;
	double mean_wei = 0;
	double mean_wei2 = 0;
	double wei = 0;
	double wei_cond = 0, mean_wei_cond = 0, mean_wei_cond2 = 0;
	for (int i = 0; i < n_sample; i++)
	{
		sample = importance_sample(delta, stat_dist, trans_mat, pwm.nrow());
		for (int j = 0; j < sample.size() - 1; j++)
		{
			sample_vec[j] = sample[j];
		}
		sample_score = compute_sample_score(pwm, sample_vec, sample[seq_len], theta);
		mean_sample += sample_score[loglik_type];
		wei = norm_const / sample_score[3];
		wei_cond = norm_const / sample_score[4] / motif_len;
		mean_wei += wei;
		mean_wei2 += wei * wei;
		mean_wei_cond += wei_cond;
		mean_wei_cond2 += wei_cond * wei_cond;
		mean_adj_score += log(sample_score[3]);
		for (int j = 0; j < scores.size(); j++)
		{
			if (scores(j) <= sample_score[loglik_type])
			{
				p_values(j, 0) += wei;
				p_values(j, 1) += wei * wei;
				p_values(j, 4) += wei_cond;
				p_values(j, 5) += wei_cond * wei_cond;
			}
		}
	}

	mean_wei /= n_sample;
	mean_wei2 /= n_sample;
	mean_wei_cond /= n_sample;
	mean_wei_cond2 /= n_sample;
	double var_wei = mean_wei2 - mean_wei * mean_wei;
	double var_wei_cond = mean_wei_cond2 - mean_wei_cond * mean_wei_cond;
	for (int j = 0; j < scores.size(); j++)
	{
		// p values
		p_values(j, 0) /= n_sample;
		p_values(j, 1) /= n_sample;
		double cov = p_values(j, 1) - mean_wei * p_values(j, 0);
		p_values(j, 1) -= p_values(j, 0) * p_values(j, 0);
		p_values(j, 2) = p_values(j, 0) / mean_wei;
		double grad1 = 1 / mean_wei;
		double grad2 = -p_values(j, 0) * grad1 * grad1;
		p_values(j, 3) = grad1 * grad1 * p_values(j, 1) + grad2 * grad2 * var_wei + 2 * grad1 * grad2 * cov;
		// weights and the weight * indicator are the same; discard the estimate
		if (var_wei == p_values(j, 1))
		{
			p_values(j, 3) = n_sample - 1;
		}
		p_values(j, 1) /= n_sample - 1;
		p_values(j, 3) /= n_sample - 1;
		// conditional p values
		p_values(j, 4) /= n_sample;
		p_values(j, 5) /= n_sample;
		cov = p_values(j, 5) - mean_wei_cond * p_values(j, 4);
		p_values(j, 5) -= p_values(j, 4) * p_values(j, 4);
		p_values(j, 6) = p_values(j, 4) / mean_wei_cond;
		grad1 = 1 / mean_wei_cond;
		grad2 = -p_values(j, 4) * grad1 * grad1;
		p_values(j, 7) = grad1 * grad1 * p_values(j, 5) + grad2 * grad2 * var_wei_cond + 2 * grad1 * grad2 * cov;
		// weights and the weight * indicator are the same; discard the estimate
		if (var_wei_cond == p_values(j, 5))
		{
			p_values(j, 7) = n_sample - 1;
		}
		p_values(j, 5) /= n_sample - 1;
		p_values(j, 7) /= n_sample - 1;
	}
	return (p_values);
}

double func_delta(NumericMatrix pwm, NumericVector stat_dist, NumericMatrix trans_mat, double theta, int seq_len)
{
	int motif_len = pwm.nrow();

	NumericMatrix delta = gen_utility_matrix(pwm, trans_mat, seq_len, theta);

	double cst = 0;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j <= seq_len - motif_len; j++)
		{
			cst += stat_dist[i] * delta(i, j);
		}
	}

	return (cst);
}

/*
Find the tilting paramter for the importance sampling distribution, using Equation (4.5).
*/
double find_theta(NumericMatrix pwm, NumericVector stat_dist, NumericMatrix trans_mat, double score, int seq_len)
{
	double theta = 0;
	double low_delta = log(func_delta(pwm, stat_dist, trans_mat, theta - 0.005, seq_len));
	double upp_delta = log(func_delta(pwm, stat_dist, trans_mat, theta + 0.005, seq_len));
	if (upp_delta - low_delta < score * 0.01)
	{
		while (upp_delta - low_delta < score * 0.01 && theta < 1)
		{
			theta += 0.01;
			low_delta = upp_delta;
			upp_delta = log(func_delta(pwm, stat_dist, trans_mat, theta + 0.005, seq_len));
		}
	}
	else
	{
		while (upp_delta - low_delta > score * 0.01 && theta > -1)
		{
			theta -= 0.01;
			upp_delta = low_delta;
			low_delta = log(func_delta(pwm, stat_dist, trans_mat, theta - 0.005, seq_len));
		}
	}
	return (theta);
}

IntegerVector importance_sample(NumericMatrix delta, NumericVector stat_dist, NumericMatrix trans_mat, int motif_len)
{
	int seq_len = delta.ncol();
	// compute the sampling distribution for each coordinate
	// sample a random vector
	RNGScope scope;
	NumericVector rv = runif(seq_len + 1);
	// note: the last digit is for sampling the start position
	// sampling the starting position of the motif
	NumericVector prob_start_pos = gen_prob_start_pos(delta, motif_len, stat_dist);

	int start_pos = sample_discrete(rv[seq_len], prob_start_pos);
	// the subsequence of length motif_len starting from start_pos follows the importance sampling distribution
	// the rest of the subsequence follows the prior distribution
	IntegerVector sample_vec(seq_len + 1);
	sample_vec[seq_len] = start_pos;

	for (int i = 0; i < seq_len; i++)
	{
		double cond_prob[4];
		for (int j = 0; j < 4; j++)
		{
			if (i == 0)
			{
				cond_prob[j] = stat_dist[j];
			}
			else
			{
				cond_prob[j] = trans_mat(sample_vec[i - 1], j);
			}
			if (i - start_pos < motif_len)
			{
				cond_prob[j] *= delta(j, seq_len - motif_len - start_pos + i);
			}
			if (j > 0)
			{
				cond_prob[j] += cond_prob[j - 1];
			}
		}
		rv[i] *= cond_prob[3];
		sample_vec[i] = 0;
		while (sample_vec[i] < 3 && rv[i] > cond_prob[sample_vec[i]])
		{
			sample_vec[i]++;
		}
	}
	return (sample_vec);
}

NumericVector compute_sample_score(NumericMatrix pwm, IntegerVector sample_vec, int start_pos, double theta)
{
	// compute the maximum score
	SequenceScores seq_scores = comp_seq_scores(pwm, sample_vec);
	// compute the weight = prior density / importance sampling density
	// NOTE: must use the score based on the true start_pos to compute the weight
	double adj_score = 0;
	for (int s = 0; s + pwm.nrow() - 1 < sample_vec.size(); s++)
	{
		adj_score += exp(theta * pwm_log_prob(pwm, sample_vec, s));
	}
	// return value
	NumericVector ret(5);
	ret[0] = seq_scores.max_log_lik;
	ret[1] = seq_scores.mean_log_lik;
	ret[2] = seq_scores.median_log_lik;
	ret[3] = adj_score;
	ret[4] = exp(theta * pwm_log_prob(pwm, sample_vec, start_pos));
	return (ret);
}

double find_percentile(NumericVector scores, double p)
{
	// compute the 1% quantile among the scores
	int n_top = scores.size() * p + 1;
	// heap stores the smalles 1% of all scores
	double heap[n_top];
	// initialize the heap
	for (int i = 0; i < n_top; i++)
	{
		heap[i] = -1e10;
	}
	// use the heap structure to find the 1% quantile among scores
	for (int i = 0; i < scores.size(); i++)
	{
		if (heap[0] < scores(i))
			heap[0] = scores(i);
		int idx = 0;
		// sort the values in the heap
		while (1)
		{
			// no children
			if (2 * idx + 1 >= n_top)
				break;
			// only one child
			if (2 * idx + 2 == n_top)
			{
				if (heap[idx] > heap[idx * 2 + 1])
				{
					double tmp = heap[idx];
					heap[idx] = heap[idx * 2 + 1];
					heap[idx * 2 + 1] = tmp;
					idx = idx * 2 + 1;
				}
				else
				{
					break;
				}
			}
			if (2 * idx + 2 < n_top)
			{
				// find the larger between the children
				int new_idx = idx * 2 + 1;
				if (heap[new_idx] > heap[new_idx + 1])
					new_idx++;
				if (heap[idx] > heap[new_idx])
				{
					double tmp = heap[idx];
					heap[idx] = heap[new_idx];
					heap[new_idx] = tmp;
					idx = new_idx;
				}
				else
				{
					break;
				}
			}
		}
	}
	return (heap[0]);
}

/* Generate a utility matrix delta for importance sampling.
Dimension of delta: 4 * seq_len.
delta(i,seq_len-1) = pwm(motif_len-1,i) ^ theta,
delta(i, j) = [sum_i' trans_mat(i, i')*delta(i, j+1)] * pwm(j-seq_len+motif_len) ^ theta,
when seq_len-motif_len <= j < seq_len-1;
delta(i, j) = [sum_i' trans_mat(i, i')*delta(i, j+1)] * pwm(j-seq_len+motif_len) ^ theta,
when 0 <= j <= seq_len-motif_len-1.
*/
NumericMatrix gen_utility_matrix(NumericMatrix pwm, NumericMatrix trans_mat, int seq_len, double theta)
{
	NumericMatrix delta(4, seq_len);
	int motif_len = pwm.nrow();
	double tol = 1e-10;
	for (int i = 0; i < 4; i++)
	{
		delta(i, seq_len - 1) = 1;
		if (theta < 0)
		{
			delta(i, seq_len - 1) = 1 / pow(pwm(motif_len - 1, i), -theta);
		}
		else
		{
			delta(i, seq_len - 1) = pow(pwm(motif_len - 1, i), theta);
		}
	}
	// Formula (A.4)
	for (int m = seq_len - 2; m >= 0; m--)
	{
		for (int i = 0; i < 4; i++)
		{
			delta(i, m) = 0;
			for (int j = 0; j < 4; j++)
			{
				delta(i, m) += trans_mat(i, j) * delta(j, m + 1);
			}
			int pos_in_motif = m - seq_len + motif_len;
			if (pos_in_motif >= 0)
			{
				if (theta < 0)
				{
					delta(i, m) /= pow(pwm(pos_in_motif, i), -theta);
				}
				else
				{
					delta(i, m) *= pow(pwm(pos_in_motif, i), theta);
				}
			}
			if (delta(i, m) < tol)
			{
				delta(i, m) = tol;
			}
		}
	}
	return delta;
}

NumericVector gen_prob_start_pos(NumericMatrix delta, int motif_len, NumericVector stat_dist)
{
	int seq_len = delta.ncol();
	NumericVector prob_start_pos(seq_len - motif_len + 1);
	for (int i = 0; i <= seq_len - motif_len; i++)
	{
		prob_start_pos[i] = 0;
		for (int j = 0; j < 4; j++)
		{
			prob_start_pos[i] += stat_dist[j] * delta(j, seq_len - motif_len - i);
		}
	}
	return prob_start_pos;
}

RcppExport SEXP test_gen_utility_matrix(SEXP _pwm, SEXP _trans_mat, SEXP _seq_len, SEXP _theta)
{
	NumericMatrix pwm(_pwm);
	NumericMatrix trans_mat(_trans_mat);
	int seq_len = as<int>(_seq_len);
	double theta = as<double>(_theta);
	return wrap(gen_utility_matrix(pwm, trans_mat, seq_len, theta));
}

RcppExport SEXP test_gen_prob_start_pos(SEXP _delta, SEXP _motif_len, SEXP _stat_dist)
{
	NumericMatrix delta(_delta);
	int motif_len = as<int>(_motif_len);
	NumericVector stat_dist(_stat_dist);
	return wrap(gen_prob_start_pos(delta, motif_len, stat_dist));
}

RcppExport SEXP test_find_percentile(SEXP _scores, SEXP _p)
{
	NumericVector scores(_scores);
	double p = as<double>(_p);
	double ret = find_percentile(scores, p);
	return (wrap(ret));
}

RcppExport SEXP compute_p_values(
	SEXP _pwm,
	SEXP _stat_dist,
	SEXP _trans_mat,
	SEXP _scores,
	SEXP _theta,
	SEXP _n_sample,
	SEXP _seq_len,
	SEXP _loglik_type)
{
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	NumericVector scores(_scores);
	double theta = as<double>(_theta);
	int n_sample = as<int>(_n_sample);
	int seq_len = as<int>(_seq_len);
	LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));

	return p_value(pwm, stat_dist, trans_mat, scores, theta, n_sample, seq_len, loglik_type);
}

RcppExport SEXP test_find_theta(SEXP _pwm, SEXP _stat_dist, SEXP _trans_mat, SEXP _score, SEXP _seq_len)
{
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	double score = as<double>(_score);
	int seq_len = as<int>(_seq_len);

	double ret = find_theta(pwm, stat_dist, trans_mat, score, seq_len);
	return (wrap(ret));
}

RcppExport SEXP test_func_delta(SEXP _pwm, SEXP _stat_dist, SEXP _trans_mat, SEXP _theta, SEXP _seq_len)
{
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	double theta = as<double>(_theta);
	int seq_len = as<int>(_seq_len);

	double ret = func_delta(pwm, stat_dist, trans_mat, theta, seq_len);
	return (wrap(ret));
}

RcppExport SEXP test_importance_sample(SEXP _delta, SEXP _stat_dist, SEXP _trans_mat, SEXP _motif_len)
{
	NumericMatrix delta(_delta);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	int motif_len = as<int>(_motif_len);
	return (wrap(importance_sample(delta, stat_dist, trans_mat, motif_len)));
}

RcppExport SEXP test_compute_sample_score(SEXP _pwm, SEXP _sample_vec, SEXP _start_pos, SEXP _theta)
{
	NumericMatrix pwm(_pwm);
	IntegerVector sample_vec(_sample_vec);
	int start_pos = as<int>(_start_pos);
	double theta = as<double>(_theta);
	return (wrap(compute_sample_score(pwm, sample_vec, start_pos, theta)));
}
