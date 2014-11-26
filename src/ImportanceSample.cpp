#include "MotifScore.h"

/*
Reference: Hock Peng Chan et al.(2010). Importance sampling of word patterns in DNA and protein sequences. Journal of computational biology, 17(12).
*/

/*
Compute the probability that a random sequence can get a score higher than 'score'.
@arg pwm The position weight matrix, with 4 columns corresponding to A, C, G, T.
@arg stat_dist A vector of length 4 with stationary distributions of A, C, G, T.
@arg trans_mat A 4 x 4 transition matrix.
@arg scores A matrix with 2 columns, with each column corresponding to one allele.
@arg p The upper percentile of the scores which is used as the mean of the importance sampling distribution.
@return A matrix with 3 columns. The first two columns are the p-values for the log-likelihood scores of each allele. The third column are the p-values for the likelihood ratios.
*/
NumericMatrix p_value(NumericMatrix pwm, NumericVector stat_dist, NumericMatrix trans_mat, NumericVector scores, double score_percentile) {
	// double score_percentile = find_percentile(scores, p);
	// printf("percentile:%3.3f\n", score_percentile);
	// find the tilting parameter
	double theta = find_theta(pwm, stat_dist, trans_mat, score_percentile);
	printf("theta:%3.3f\n", theta);
	NumericMatrix p_values(scores.size(), 4);

	double tol = 1e-10;
	int motif_len = pwm.nrow();
	NumericMatrix delta(4, motif_len * 2 - 1);
	IntegerVector sample_vec(2 * motif_len - 1);
	IntegerVector sample(2 * motif_len);
	NumericVector sample_score(5);

	for(int i = 0; i < 4; i ++)
		for(int m = 0; m < motif_len; m ++)
			if(pwm(m, i) < tol)
				pwm(m, i) = tol;
	for(int i = 0; i < 4; i ++) {
		delta(i, 2 * motif_len - 2) = 1;
		if(theta < 0) {
			delta(i, 2 * motif_len - 2) = 1 / pow(pwm(motif_len - 1, i), - theta);
		} else {
			delta(i, 2 * motif_len - 2) = pow(pwm(motif_len - 1, i), theta);
		}
		//printf("delta(%d,%d)=%3.10f\n", i, motif_len - 1, delta(i, motif_len - 1));
	}
	// Formula (A.4)
	for(int m = 2 * motif_len - 3; m >= 0; m --) {
		for(int i = 0; i < 4; i ++) {
			delta(i, m) = 0;
			for(int j = 0; j < 4; j ++) {
				delta(i, m) += trans_mat(i, j) * delta(j, m + 1);
			}
			if(m >= motif_len - 1) {
				if(theta < 0) {
					delta(i, m) /= pow(pwm(m - motif_len + 1, i), - theta);
				} else {
					delta(i, m) *= pow(pwm(m - motif_len + 1, i), theta);
				}
			}
			if(delta(i, m) < tol) {
				delta(i, m) = tol;
			}
			//printf("delta(%d,%d)=%3.10f\n", i, m, delta(i, m));
		}
	}

	double norm_const = 0;
	for(int i = 0; i < 4; i ++) {
		for(int j = 0; j < motif_len; j ++) {
			norm_const += stat_dist[i] * delta(i, j);
		}
	}
	//	printf("Constant value : %3.10f\n", norm_const);

	for(int i = 0; i < p_values.nrow(); i ++)
		for(int j = 0; j < 4; j ++)
			p_values(i, j) = 0;
	
	int n_sample = 1e4;
	double mean_sample = 0;
	double mean_adj_score = 0;
	double mean_wei = 0;
	double mean_wei2 = 0;
	double wei = 0;
	for(int i = 0; i < n_sample; i ++) {
		sample = importance_sample(delta, stat_dist, trans_mat, pwm, theta);
		for(int j = 0; j < sample.size() - 1; j ++) {
			sample_vec[j] = sample[j];
		}
		sample_score = compute_sample_score(pwm, sample_vec, sample[2 * motif_len - 1], theta);
		mean_sample += sample_score[0];
		wei = norm_const / sample_score[1];
		mean_wei += wei;
		mean_wei2 += wei * wei;
		mean_adj_score += log(sample_score[1]);
		for(int j = 0; j < scores.size(); j ++) {
			if(scores(j) <= sample_score[0]) {
				p_values(j, 0) += wei;
				p_values(j, 1) += wei * wei;
			}
		}
	}
	printf("Mean sample : %lf \t adj_score : %lf \t weight : %lf \n", mean_sample / n_sample, mean_adj_score / n_sample, mean_wei / n_sample);
	mean_wei /= n_sample;
	mean_wei2 /= n_sample;
	double var_wei = mean_wei2 - mean_wei * mean_wei;
	for(int j = 0; j < scores.size(); j ++) {
		p_values(j, 0) /= n_sample;
		p_values(j, 1) /= n_sample;
		double cov = p_values(j, 1) - mean_wei * p_values(j, 0);
		p_values(j, 1) -= p_values(j, 0) * p_values(j, 0);
		p_values(j, 2) = p_values(j, 0) / mean_wei;
		double grad1 = 1 / mean_wei;
		double grad2 = - p_values(j, 0) * grad1 * grad1;
		p_values(j, 3) = grad1 * grad1 * p_values(j, 1) + grad2 * grad2 * var_wei + 2 * grad1 * grad2 * cov;
		if(var_wei == p_values(j, 1)) {
			p_values(j, 3) = n_sample - 1;
		}
		p_values(j, 1) /= n_sample - 1;
		p_values(j, 3) /= n_sample - 1;
	}
	return(p_values);
}

double func_delta(NumericMatrix pwm, NumericVector stat_dist, NumericMatrix trans_mat, double theta) {
	int motif_len = pwm.nrow();
	double tol = 1e-10;

	NumericMatrix delta(4, motif_len * 2 - 1);

	for(int i = 0; i < 4; i ++)
		for(int m = 0; m < motif_len; m ++)
			if(pwm(m, i) < tol)
				pwm(m, i) = tol;
	for(int i = 0; i < 4; i ++) {
		delta(i, 2 * motif_len - 2) = 1;
		if(theta < 0) {
			delta(i, 2 * motif_len - 2) = 1 / pow(pwm(motif_len - 1, i), - theta);
		} else {
			delta(i, 2 * motif_len - 2) = pow(pwm(motif_len - 1, i), theta);
		}
		//printf("delta(%d,%d)=%3.10f\n", i, motif_len - 1, delta(i, motif_len - 1));
	}
	// Formula (A.4)
	for(int m = 2 * motif_len - 3; m >= 0; m --) {
		for(int i = 0; i < 4; i ++) {
			delta(i, m) = 0;
			for(int j = 0; j < 4; j ++) {
				delta(i, m) += trans_mat(i, j) * delta(j, m + 1);
			}
			if(m >= motif_len - 1) {
				if(theta < 0) {
					delta(i, m) /= pow(pwm(m - motif_len + 1, i), - theta);
				} else {
					delta(i, m) *= pow(pwm(m - motif_len + 1, i), theta);
				}
			}
			if(delta(i, m) < tol) {
				delta(i, m) = tol;
			}
			//printf("delta(%d,%d)=%3.10f\n", i, m, delta(i, m));
		}
	}

	double cst = 0;
	for(int i = 0; i < 4; i ++) {
		for(int j = 0; j < motif_len; j ++) {
			cst += stat_dist[i] * delta(i, j);
		}
	}
	
	return(cst);
}

/*
Find the tilting paramter for the importance sampling distribution, using Equation (4.5).
*/
double find_theta(NumericMatrix pwm, NumericVector stat_dist, NumericMatrix trans_mat, double score) {
	double theta = 0;
	double low_delta = log(func_delta(pwm, stat_dist, trans_mat, theta - 0.005));
	double upp_delta = log(func_delta(pwm, stat_dist, trans_mat, theta + 0.005));
	if(upp_delta - low_delta < score * 0.01) {
		while(upp_delta - low_delta < score * 0.01 && theta < 1) {
			theta += 0.01;
			low_delta = upp_delta;
			upp_delta = log(func_delta(pwm, stat_dist, trans_mat, theta + 0.005));
		}
	} else {
		while(upp_delta - low_delta > score * 0.01 && theta > -1) {
			theta -= 0.01;
			upp_delta = low_delta;
			low_delta = log(func_delta(pwm, stat_dist, trans_mat, theta - 0.005));
		}
	}
	return(theta);
}

IntegerVector importance_sample(NumericMatrix delta, NumericVector stat_dist, NumericMatrix trans_mat, NumericMatrix pwm, double theta) {
	int motif_len = pwm.nrow();
	// compute the sampling distribution for each coordinate
	// sample a random vector
	RNGScope scope;
	NumericVector rv = runif(2 * motif_len);
	// note: the last digit is for sampling the start position
	// sampling the starting position of the motif
	double prob_stat[motif_len];
	for(int i = 0; i < motif_len; i ++) {
		prob_stat[motif_len - 1 - i] = 0;
		for(int j = 0; j < 4; j ++) {
			prob_stat[motif_len - i - 1] += stat_dist[j] * delta(j, i);
		}
		if(i > 0)
			prob_stat[i] += prob_stat[i - 1];
	}
	
	rv[2 * motif_len - 1] *= prob_stat[motif_len - 1];
	int start_pos = 0;
	while(rv[2 * motif_len - 1] > prob_stat[start_pos]) {
		start_pos ++;
	}
	// the subsequence of length motif_len starting from start_pos follows the importance sampling distribution
	// the rest of the subsequence follows the prior distribution
	IntegerVector sample_vec(motif_len * 2);
	sample_vec[motif_len * 2 - 1] = start_pos;
	
	for(int i = 0; i < 2 * motif_len - 1; i ++) {
		double cond_prob[4];
		for(int j = 0; j < 4; j ++) {
			if(i == 0) {
				cond_prob[j] = stat_dist[j];
			} else {
				cond_prob[j] = trans_mat(sample_vec[i - 1], j);
			}
			if(motif_len - 1 - start_pos + i < motif_len * 2 - 1) {
				cond_prob[j] *= delta(j, motif_len - 1 - start_pos + i);
			}
			if(j > 0) {
				cond_prob[j] += cond_prob[j - 1];
			}
		}
		/*
		if(i >= start_pos && i < start_pos + motif_len) {
			double test_prob[4];
			for(int j = 0; j < 4; j ++) {
				test_prob[j] = stat_dist[j];
				test_prob[j] *= pow(pwm(i - start_pos, j), theta);
				if(j > 0)
					test_prob[j] += test_prob[j - 1];
			}
			for(int j = 0; j < 4; j ++) {
				printf("%d,%d: prob = %3.10f \t %3.10f \n", i - start_pos, j, cond_prob[j] / cond_prob[3], test_prob[j] / test_prob[3]);
			}
		}
		*/
		rv[i] *= cond_prob[3];
		sample_vec[i] = 0;
		while(sample_vec[i] < 3 && rv[i] > cond_prob[sample_vec[i]]) {
			sample_vec[i] ++;
		}
		// index from 1
		// sample_vec[i] ++;
		//		printf("%d\t", sample_vec[i]);
	}
	return(sample_vec);
}

NumericVector compute_sample_score(NumericMatrix pwm, IntegerVector sample_vec, int start_pos, double theta) {
	int seq_len = sample_vec.size();
	//	printf("\n");
	// compute the reverse strand sequence
	IntegerVector rev_sample_vec(seq_len);
	IntegerVector sample_vec_copy(seq_len);
	IntegerVector rev_sample_vec_copy(seq_len);
	for(int i = 0; i < seq_len; i ++) {
		rev_sample_vec[i] = 3 - sample_vec[seq_len - 1 - i];
		sample_vec_copy[i] = sample_vec[i];
		rev_sample_vec_copy[i] = 3 - sample_vec[seq_len - 1 - i];
	}
	// compute the maximum score
	double rnd_score = pwm_log_prob(pwm, sample_vec, find_best_match(pwm, sample_vec));
	double rnd_score_rev = pwm_log_prob(pwm, rev_sample_vec, find_best_match(pwm, rev_sample_vec));
	if(rnd_score_rev > rnd_score)
		rnd_score = rnd_score_rev;
	// SNP score
	//	double snp_score[3];
	int snp_id = 0;
	double rnd_score_copy, rnd_score_rev_copy;
	for(int j = 0; j < 4; j ++) {
		if(sample_vec[seq_len / 2] == j)
			continue;
		sample_vec_copy[seq_len / 2] = j;
		rev_sample_vec_copy[seq_len / 2] = 3 - j;
		rnd_score_copy = pwm_log_prob(pwm, sample_vec_copy, find_best_match(pwm, sample_vec_copy));
		rnd_score_rev_copy = pwm_log_prob(pwm, rev_sample_vec_copy, find_best_match(pwm, rev_sample_vec_copy));
		if(rnd_score_rev_copy > rnd_score_copy)
			rnd_score_copy = rnd_score_rev_copy;
		//		snp_score[snp_id] = rnd_score_copy - rnd_score;
		snp_id ++;
	}
	// compute the weight = prior density / importance sampling density
	// note: must use the score based on the true start_pos to compute the weight
	// this is a bug that took 2 days to fix!
	double adj_score = 0;
	for(int s = 0; s < pwm.nrow(); s ++) {
		adj_score += exp(theta * pwm_log_prob(pwm, sample_vec, s));
	}
	// return value
	NumericVector ret(2);
	ret[0] = rnd_score;
	ret[1] = adj_score;
	//	ret[2] = snp_score[0];
	//	ret[3] = snp_score[1];
	//	ret[4] = snp_score[2];
	//	printf("score:%3.3f\tweight:%3.3f\tconstant:%3.3f\n", ret[0], ret[1], log(delta(0, 3)));
	return(ret);
	
}

double find_percentile(NumericVector scores, double p) {
	// compute the 1% quantile among the scores
	int n_top = scores.size() * p + 1;
	// heap stores the smalles 1% of all scores
	double heap[n_top];
	// initialize the heap
	for(int i = 0; i < n_top; i ++) {
		heap[i] = -1e10;
	}
	// use the heap structure to find the 1% quantile among scores
	for(int i = 0; i < scores.size(); i ++) {
		if(heap[0] < scores(i))
			heap[0] = scores(i);
		int idx = 0;
		// sort the values in the heap
		while(1) {
			// no children
			if(2 * idx + 1 >= n_top)
				break;
			// only one child
			if(2 * idx + 2 == n_top) {
				if(heap[idx] > heap[idx * 2 + 1]) {
					double tmp = heap[idx];
					heap[idx] = heap[idx * 2 + 1];
					heap[idx * 2 + 1] = tmp;
					idx = idx * 2 + 1;
				} else {
					break;
				}
			}
			if(2 * idx + 2 < n_top) {
				// find the larger between the children
				int new_idx = idx * 2 + 1;
				if(heap[new_idx] > heap[new_idx + 1])
					new_idx ++;
				if(heap[idx] > heap[new_idx]) {
					double tmp = heap[idx];
					heap[idx] = heap[new_idx];
					heap[new_idx] = tmp;
					idx = new_idx;
				} else {
					break;
				}
			}
		}
	}
	return(heap[0]);
}

SEXP test_find_percentile(SEXP _scores, SEXP _p) {
	NumericVector scores(_scores);
	double p = as<double>(_p);
	double ret = find_percentile(scores, p);
	return(wrap(ret));
}

SEXP test_p_value(SEXP _pwm, SEXP _stat_dist, SEXP _trans_mat, SEXP _scores, SEXP _perc) {
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	NumericVector scores(_scores);
	double perc = as<double>(_perc);
	
	NumericMatrix p_values = p_value(pwm, stat_dist, trans_mat, scores, perc);
	return(wrap(p_values));
}

SEXP test_find_theta(SEXP _pwm, SEXP _stat_dist, SEXP _trans_mat, SEXP _score) {
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	double score = as<double>(_score);

	double ret = find_theta(pwm, stat_dist, trans_mat, score);
	return(wrap(ret));
}

SEXP test_func_delta(SEXP _pwm, SEXP _stat_dist, SEXP _trans_mat, SEXP _theta) {
	NumericMatrix pwm(_pwm);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	double theta = as<double>(_theta);
	
	double ret = func_delta(pwm, stat_dist, trans_mat, theta);
	return(wrap(ret));
}

SEXP test_importance_sample(SEXP _delta, SEXP _stat_dist, SEXP _trans_mat, SEXP _pwm, SEXP _theta) {
	NumericMatrix delta(_delta);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	NumericMatrix pwm(_pwm);
	double theta = as<double>(_theta);
	return(wrap(importance_sample(delta, stat_dist, trans_mat, pwm, theta)));
}

SEXP test_compute_sample_score(SEXP _pwm, SEXP _sample_vec, SEXP _start_pos, SEXP _theta) {
	NumericMatrix pwm(_pwm);
	IntegerVector sample_vec(_sample_vec);
	int start_pos = as<int>(_start_pos);
	double theta = as<double>(_theta);
	return(wrap(compute_sample_score(pwm, sample_vec, start_pos, theta)));
}
