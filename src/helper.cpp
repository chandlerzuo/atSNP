#include "helper.h"

int sample_discrete(double rnd, NumericVector pdf) {
    assert(rnd>=0 && rnd <=1);
    double norm_const = 0;
    for(int i = 0; i < pdf.size(); ++ i) {
        norm_const += pdf[i];
    }
    rnd *= norm_const;
    int pos = 0;
    double cum_sum = pdf[0];
    while(cum_sum < rnd) {
        assert(pos < pdf.size() - 1);
        cum_sum += pdf[++pos];
    }
    return pos;
}

void rowwise_l1_normalize(NumericMatrix & mat, double epsilon) {
    for(int i = 0; i < mat.nrow(); ++ i) {
        double row_sum = 0;
        for(int j = 0; j < mat.ncol(); ++ j) {
            if(mat(i, j) < epsilon) {
                mat(i, j) = epsilon;
            }
            row_sum += mat(i, j);
        }
        for(int j = 0; j < mat.ncol(); ++ j) {
            mat(i, j) /= row_sum;
        }
    }
}


/*
@arg scores A list of scores.
@arg weights A list of weights.
@arg sample_score A matrix with Monte Carlo sampled scores. The number 
of rows is the sample size. The number of columns is the number of 
Monte Carlo examples for each example. For Indel, ncol=1. For SNP, 
ncol=3.
@return A matrix of 4 columns. Column 0: p-value based on moment estimates.
Column 1: variance of Column 0. Column 2: p-value based on ratio estimates.
Column 3: variance of Column 2.
*/
NumericMatrix comp_empirical_p_values(
    NumericVector scores,
    NumericVector weights,
    NumericMatrix sample_score
) {
	double wei_sum = 0;
	double wei2_sum = 0;
	int sample_size = sample_score.nrow();
    assert(weights.size() == sample_size);
	int n_scores = scores.size();
	NumericMatrix p_values(n_scores, 4);
	NumericMatrix moments(n_scores, 4);
	for(int j = 0; j < n_scores; j ++) {
		for(int k = 0; k < 2; k ++) {
			moments(j, k) = 0;
		}
	}
	for(int i = 0; i < sample_size; i ++) {
		wei_sum += weights[i];
		wei2_sum += weights[i] * weights[i];
		for(int j = 0; j < n_scores; j ++) {
			for(int k = 0; k < sample_score.ncol(); k ++) {
				// SNP changes binding affinity
				if(sample_score(i, k) >= scores(j) || sample_score(i, k) <= -scores(j)) {
					moments(j, 0) += weights[i];
					moments(j, 1) += weights[i] * weights[i];
				}
			}
		}
	}
	wei2_sum /= sample_size;
	wei_sum /= sample_size;
	double var2 = wei2_sum - wei_sum * wei_sum;
	double grad1 = 1 / wei_sum;
	double grad2 = 0, var1 = 0, cov = 0;
	for(int j = 0; j < n_scores; j ++) {
		p_values(j, 0) = moments(j, 0) / 3 / sample_size;
		p_values(j, 1) = moments(j, 1) / 3 / sample_size - p_values(j, 0) * p_values(j, 0);
		p_values(j, 2) = p_values(j, 0) / wei_sum;
		grad2 = - p_values(j, 0) * grad1 * grad1;
		var1 = p_values(j, 1);
		cov = moments(j, 1) / 3 / sample_size - p_values(j, 0) * wei_sum;
		p_values(j, 1) /= 3 * sample_size - 1;
		if(p_values(j, 0) != wei_sum) {
			p_values(j, 3) = grad1 * grad1 * var1 + grad2 * grad2 * var2 + 2 * grad1 * grad2 * cov;
			p_values(j, 3) /= 3 * sample_size - 1;
		} else {
			p_values(j, 3) = 1;
		}
	}
	return(p_values);
}
