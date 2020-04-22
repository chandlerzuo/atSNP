#include "MotifScore.h"
#include "helper.h"

/*
Reference: Hock Peng Chan et al.(2010). Importance sampling of word patterns in DNA and protein sequences. Journal of computational biology, 17(12).
*/

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
Rcpp::List p_value_change(
  NumericMatrix pwm,
  NumericMatrix wei_mat,
  NumericMatrix adj_mat,
  NumericVector stat_dist,
  NumericMatrix trans_mat,
  NumericVector scores,
  NumericVector pval_ratio,
  double score_percentile,
  int n_sample,
  LoglikType loglik_type
) {
	// find the tilting parameter
	double theta = find_theta_change(wei_mat, adj_mat, score_percentile);
	NumericMatrix p_values(scores.size(), 4);
	NumericVector sample_score(5);
	int start_pos;

	double tol = 1e-10;
	int motif_len = pwm.nrow();
	IntegerVector sample(2 * motif_len);
	IntegerVector sample_vec(2 * motif_len - 1);

	for(int i = 0; i < 4; i ++)
		for(int m = 0; m < motif_len; m ++) {
			if(pwm(m, i) < tol)
				pwm(m, i) = tol;
			if(wei_mat(m, i) < tol)
				wei_mat(m, i) = tol;
		}
	

	double norm_const = func_delta_change(wei_mat, adj_mat, theta);

	for(int i = 0; i < p_values.nrow(); i ++) {
		for(int j = 0; j < 4; j ++) {
			p_values(i, j) = 0;
		}
	}	
	double wei = 0;
	double mean_score = 0;
	NumericVector weights(n_sample);
	double score[n_sample][4];
	NumericMatrix score_diff_sam(n_sample, 3);
	for(int i = 0; i < n_sample; i ++) {
		sample = importance_sample_change(
		  adj_mat, stat_dist, trans_mat, wei_mat, theta
	  );
		for(int j = 0; j < motif_len * 2 - 1; j ++) {
			sample_vec[j] = sample[j];
		}
		start_pos = sample[2 * motif_len - 1];
		sample_score = compute_sample_score_change(
		  pwm,
		  wei_mat,
		  adj_mat,
		  sample_vec,
		  stat_dist,
		  trans_mat,
		  start_pos,
		  theta,
		  loglik_type
		);
		wei = norm_const / sample_score[0];
		mean_score += sample_score[4];
		// copy the weights and the scores for each allele
		score[i][0] = sample_score[5];
		weights[i] = wei;
		for(int j = 0; j < 3; j ++) {
			score[i][j + 1] = sample_score[5] - sample_score[j + 1];
			if(sample_score[j + 1] >= 0) {
				score_diff_sam(i, j) = sample_score[j + 1];
			} else {
				score_diff_sam(i, j) = - sample_score[j + 1];
			}
		}
	}

	NumericMatrix pval_loglik(scores.size(), 3);
	pval_loglik = comp_empirical_p_values(scores, weights, score_diff_sam);
	
	// compute the sample log ranks
	double pval_sam[4];
	NumericMatrix pval_ratio_sam(n_sample, 3);
	for(int i = 0; i < n_sample; i ++) {
		for(int j = 0; j < 4; j ++) {
			pval_sam[j] = 0;
		}
		for(int i1 = 0; i1 < n_sample; i1 ++) {
			for(int j = 0; j < 4; j ++) {
				if(score[i1][0] >= score[i][j]) {
					pval_sam[j] += weights[i1];
				}
			}
		}
		for(int j = 0; j < 4; j ++) {
			if(pval_sam[j] < tol) {
				pval_sam[j] = tol;
			}
		}
		for(int j = 0; j < 3; j ++) {
			pval_ratio_sam(i, j) = log(pval_sam[0]) - log(pval_sam[j + 1]);
			if(pval_ratio_sam(i, j) < 0) {
				pval_ratio_sam(i, j) = - pval_ratio_sam(i, j);
			}
		}
	}

	NumericMatrix pval_rank(scores.size(), 4);
	pval_rank = comp_empirical_p_values(pval_ratio, weights, pval_ratio_sam);

	Rcpp::List ret = Rcpp::List::create(
	  Rcpp::Named("score") = pval_loglik,
	  Rcpp::Named("rank") = pval_rank
	);
	return(ret);
}

double func_delta_change(NumericMatrix wei_mat, NumericMatrix adj_mat, double theta) {
	int motif_len = wei_mat.nrow();
	NumericMatrix delta(4, motif_len);
	
	double adj_sum[motif_len];
	for(int i = 0; i < motif_len; i ++) {
		adj_sum[i] = 0;
		for(int j = 0; j < 4; j ++) {
			adj_sum[i] += adj_mat(i, j);
		}
	}

	double norm_const = 0;
	double tmp = 0;
	for(int i = 0; i < motif_len; i ++) {
		tmp = 0;
		for(int j = 0; j < 4; j ++) {
			tmp += exp(log(wei_mat(i, j)) * theta);
		}
		// two stationary sequences, one with length i, one with length motif_len - i - 1
		norm_const += tmp / adj_sum[i];
	}

	for(int i = 0; i < motif_len; i ++) {
		norm_const *= adj_sum[i];
	}

	return(norm_const);
}

/*
  Find the tilting paramter for the importance sampling distribution, using Equation (4.5).
*/
double find_theta_change(NumericMatrix wei_mat, NumericMatrix adj_mat, double score) {
	double theta = 0;
	double tol = 0.01;
	double low_delta = log(func_delta_change(wei_mat, adj_mat, theta - tol / 2));
	double upp_delta = log(func_delta_change(wei_mat, adj_mat, theta + tol / 2));
	if(upp_delta - low_delta < score * tol) {
		while(upp_delta - low_delta < score * tol && theta < 1) {
			theta += tol;
			low_delta = upp_delta;
			upp_delta = log(func_delta_change(wei_mat, adj_mat, theta + tol / 2));
			if(upp_delta - low_delta >= score * tol && tol > 1e-4) {
				theta -= tol;
				tol /= 10;
				low_delta = log(func_delta_change(wei_mat, adj_mat, theta - tol / 2));
				upp_delta = log(func_delta_change(wei_mat, adj_mat, theta + tol / 2));;
			}
		}
	} else {
		while(upp_delta - low_delta > score * tol && theta > -1) {
			theta -= tol;
			upp_delta = low_delta;
			low_delta = log(func_delta_change(wei_mat, adj_mat, theta - tol / 2));
			if(upp_delta - low_delta <= score * tol && tol > 1e-4) {
				theta += tol;
				tol /= 10;
				low_delta = log(func_delta_change(wei_mat, adj_mat, theta - tol / 2));
				upp_delta = log(func_delta_change(wei_mat, adj_mat, theta + tol / 2));;
			}
		}
	}
	return(theta);
}

IntegerVector importance_sample_change(NumericMatrix adj_mat, NumericVector stat_dist, NumericMatrix trans_mat, NumericMatrix wei_mat, double theta) {
	int motif_len = wei_mat.nrow();
	// compute the sampling distribution for each coordinate
	// sample a random vector
	RNGScope scope;
	NumericVector rv = runif(2 * motif_len);
	// note: the last digit is for sampling the start position
	// 1. sample the starting position
	double prob_start[motif_len];
	double tmp;
	for(int i = 0; i < motif_len; i ++) {
		prob_start[i] = 0;
		tmp = 0;
		for(int j = 0; j < 4; j ++) {
			prob_start[i] += exp(log(wei_mat(motif_len - 1 - i, j)) * theta);
			tmp += adj_mat(motif_len - 1 - i, j);
		}
		prob_start[i] /= tmp;
		if(i > 0)
			prob_start[i] += prob_start[i - 1];
	}
	double norm_const = prob_start[motif_len - 1];
	rv[2 * motif_len - 1] *= norm_const;
	int start_pos = 0;
	while(start_pos < motif_len - 1 && rv[2 * motif_len - 1] > prob_start[start_pos]) {
		start_pos ++;
	}
	// the rest of the subsequence follows the prior distribution
	if(start_pos == motif_len)
		start_pos = motif_len - 1;
	// 2. sample the actual vector
	IntegerVector sample_vec(motif_len * 2);
	sample_vec[motif_len * 2 - 1] = start_pos;
	for(int i = 0; i < 2 * motif_len - 1; i ++) {
		if(i < start_pos || i > start_pos + motif_len - 1) {
			double cond_prob[4];
			for(int j = 0; j < 4; j ++) {
				if(i == 0) {
					cond_prob[j] = stat_dist[j];
				} else {
					cond_prob[j] = trans_mat(sample_vec[i - 1], j);
				}
				if(j > 0) {
					cond_prob[j] += cond_prob[j - 1];
				}
			}
			rv[i] *= cond_prob[3];
			sample_vec[i] = 0;
			while(sample_vec[i] < 3 && rv[i] > cond_prob[sample_vec[i]]) {
				sample_vec[i] ++;
			}
		}
	}
	for(int i = start_pos; i < start_pos + motif_len; i ++) {
		double cond_prob[4];
		for(int j = 0; j < 4; j ++) {
			if(i == motif_len - 1) {
				cond_prob[j] = exp(log(wei_mat(i - start_pos, j)) * theta);
			} else {
				cond_prob[j] = adj_mat(i - start_pos, j);
			}
			if(j > 0) {
				cond_prob[j] += cond_prob[j - 1];
			}
		} 
		rv[i] *= cond_prob[3];
		sample_vec[i] = 0;
		while(sample_vec[i] < 3 && rv[i] > cond_prob[sample_vec[i]]) {
			sample_vec[i] ++;
		}
	}
	return(sample_vec);
}

NumericVector compute_sample_score_change(
  NumericMatrix pwm,
  NumericMatrix wei_mat,
  NumericMatrix adj_mat,
  IntegerVector sample_vec,
  NumericVector stat_dist,
  NumericMatrix trans_mat,
  int start_pos,
  double theta,
  LoglikType loglik_type
) {
	// compute the reverse strand sequence
	int motif_len = pwm.nrow();
	IntegerVector sample_vec_copy(motif_len * 2 - 1);
	for(int i = 0; i < motif_len * 2 - 1; i ++) {
		sample_vec_copy[i] = sample_vec[i];
	}
	// compute the maximum score
	double rnd_score = 0;
	SequenceScores seq_scores = comp_seq_scores(pwm, sample_vec);
	switch(loglik_type){
	case LoglikType::mean:
	  rnd_score = seq_scores.mean_log_lik;
	  break;
	case LoglikType::median:
	  rnd_score = seq_scores.median_log_lik;
	  break;
	default:
	  rnd_score = seq_scores.max_log_lik;
	  break;
	}
	// SNP score
	double snp_score[3];
	int snp_id = 0;
	double rnd_score_copy;
	for(int j = 0; j < 4; j ++) {
		if(sample_vec[motif_len - 1] == j)
			continue;
		sample_vec_copy[motif_len - 1] = j;
		SequenceScores seq_scores_for_copy = comp_seq_scores(pwm, sample_vec_copy);
		switch(loglik_type){
		case LoglikType::mean:
		  rnd_score_copy = seq_scores_for_copy.mean_log_lik;
		  break;
		case LoglikType::median:
		  rnd_score_copy = seq_scores_for_copy.median_log_lik;
		  break;
		default:
		  rnd_score_copy = seq_scores_for_copy.max_log_lik;
  		break;
		}
		snp_score[snp_id] = rnd_score - rnd_score_copy;
		snp_id ++;
	}
	if(snp_id != 3){
		throw std::length_error("snp_score has a wrong size.");
	};

	// compute the weight = prior density / importance sampling density
	// note: must use the score based on the true start_pos to compute the weight
	// this is a bug that took 2 days to fix!
	double adj_score = 0;
	double adj_s = 0;
	for(int s = 0; s < motif_len; s ++) {
		adj_s = 0;
		for(int i = 0; i < motif_len; i ++) {
			if(i == motif_len - 1 - s) {
				adj_s += theta * log(wei_mat(motif_len - 1 - s, sample_vec[motif_len - 1]));
			} else {
				adj_s += log(adj_mat(i, sample_vec[s + i]));
			}
			if(s + i == 0) {
				adj_s -= log(stat_dist[sample_vec[s + i]]);
			} else {
				adj_s -= log(trans_mat(sample_vec[s + i - 1], sample_vec[s + i]));
			}
		}
		if(s < motif_len - 1) {
			adj_s += log(stat_dist[sample_vec[s + motif_len]]) -
				log(trans_mat(sample_vec[s + motif_len - 1], sample_vec[s + motif_len]));
		}
		adj_score += exp(adj_s);
	}
	// return value
	NumericVector ret(6);
	ret[0] = adj_score;
	ret[1] = snp_score[0];
	ret[2] = snp_score[1];
	ret[3] = snp_score[2];
	ret[4] = log(wei_mat(motif_len - 1 - start_pos, sample_vec[motif_len - 1]));
	ret[5] = rnd_score;
	return(ret);
}

double find_percentile_change(NumericVector scores, double p) {
	// compute the 1% quantile among the scores
	int n_top = scores.size() * p + 1;
	// heap stores the smalles 1% of all scores
	double heap[n_top];
	// initialize the heap
	for(int i = 0; i < n_top; i ++) {
		heap[i] = -1e10;
	}
	double score_change = 0;
	// use the heap structure to find the 1% quantile among scores
	for(int i = 0; i < scores.size(); i ++) {
		score_change = scores(i);
		if(score_change < 0)
			score_change = - score_change;
		if(heap[0] < score_change)
			heap[0] = score_change;
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

RcppExport SEXP test_find_percentile_change(SEXP _scores, SEXP _perc) {
	NumericVector scores(_scores);
	double perc = as<double>(_perc);
	double ret = find_percentile_change(scores, perc);
	return(wrap(ret));
}

SEXP compute_p_value_change(
  SEXP _pwm,
  SEXP _wei_mat,
  SEXP _adj_mat,
  SEXP _stat_dist,
  SEXP _trans_mat,
  SEXP _scores,
  SEXP _pval_ratio,
  SEXP _perc,
  SEXP _n_sample,
  SEXP _loglik_type
) {
	NumericMatrix pwm(_pwm);
	NumericMatrix wei_mat(_wei_mat);
	NumericMatrix adj_mat(_adj_mat);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	NumericVector scores(_scores);
	NumericVector pval_ratio(_pval_ratio);
	double perc = as<double>(_perc);
	int n_sample = as<int>(_n_sample);
	LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));

	Rcpp::List p_values = p_value_change(
	  pwm,
	  wei_mat,
	  adj_mat,
	  stat_dist,
	  trans_mat,
	  scores,
	  pval_ratio,
	  perc,
	  n_sample,
	  loglik_type
	 );
	return(wrap(p_values));
}

SEXP test_find_theta_change(SEXP _wei_mat, SEXP _adj_mat, SEXP _score) {
	NumericMatrix wei_mat(_wei_mat);
	NumericMatrix adj_mat(_adj_mat);
	double score = as<double>(_score);

	double ret = find_theta_change(wei_mat, adj_mat, score);
	return(wrap(ret));
}

SEXP test_func_delta_change(SEXP _wei_mat, SEXP _adj_mat, SEXP _theta) {
	NumericMatrix wei_mat(_wei_mat);
	NumericMatrix adj_mat(_adj_mat);
	double theta = as<double>(_theta);
	
	double ret = func_delta_change(wei_mat, adj_mat, theta);
	return(wrap(ret));
}

SEXP test_importance_sample_change(SEXP _adj_mat, SEXP _stat_dist, SEXP _trans_mat, SEXP _wei_mat, SEXP _theta) {
	NumericMatrix adj_mat(_adj_mat);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	NumericMatrix wei_mat(_wei_mat);
	double theta = as<double>(_theta);
	return(wrap(importance_sample_change(adj_mat, stat_dist, trans_mat, wei_mat, theta)));
}

SEXP test_compute_sample_score_change(
  SEXP _pwm,
  SEXP _wei_mat,
  SEXP _adj_mat,
  SEXP _sample_vec,
  SEXP _stat_dist,
  SEXP _trans_mat,
  SEXP _start_pos,
  SEXP _theta,
  SEXP _loglik_type
) {
	NumericMatrix pwm(_pwm);
	NumericMatrix wei_mat(_wei_mat);
	NumericMatrix adj_mat(_adj_mat);
	NumericVector stat_dist(_stat_dist);
	NumericMatrix trans_mat(_trans_mat);
	IntegerVector sample_vec(_sample_vec);
	int start_pos = as<int>(_start_pos);
	double theta = as<double>(_theta);
	LoglikType loglik_type = static_cast<LoglikType>(as<int>(_loglik_type));
	return(
	  wrap(
	    compute_sample_score_change(
	      pwm,
	      wei_mat,
	      adj_mat,
  	    sample_vec,
	      stat_dist,
	      trans_mat,
	      start_pos,
	      theta,
	      loglik_type
      )
	  )
	);
}
