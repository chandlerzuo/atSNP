#include "helper.h"

int sample_discrete(double rnd, NumericVector pdf)
{
	if (rnd < 0 || rnd > 1)
	{
		throw std::invalid_argument("Random number has to be between 0 and 1.");
	};
	double norm_const = 0;
	for (int i = 0; i < pdf.size(); ++i)
	{
		norm_const += pdf[i];
	}
	rnd *= norm_const;
	int pos = 0;
	double cum_sum = pdf[0];
	while (cum_sum < rnd)
	{
		if (pos >= pdf.size() - 1)
		{
			throw std::range_error("Index out of range.");
		};
		cum_sum += pdf[++pos];
	}
	return pos;
}

void rowwise_l1_normalize(NumericMatrix &mat, double epsilon)
{
	for (int i = 0; i < mat.nrow(); ++i)
	{
		double row_sum = 0;
		for (int j = 0; j < mat.ncol(); ++j)
		{
			if (mat(i, j) < epsilon)
			{
				mat(i, j) = epsilon;
			}
			row_sum += mat(i, j);
		}
		for (int j = 0; j < mat.ncol(); ++j)
		{
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
@arg test_type Test type. It defines the rejection region as where
gte/lte: Monte Carlo example statistics are gte/lte than scores;
abs: the abs of Monte Carlo statistics are gte than scores;
two_sided: 2 x the smaller of rejection regions by gte or lte.
@return A matrix of 4 columns. Column 0: p-value based on moment estimates.
Column 1: variance of Column 0. Column 2: p-value based on ratio estimates.
Column 3: variance of Column 2.
*/
NumericMatrix comp_empirical_p_values(
	NumericVector scores,
	NumericVector weights,
	NumericMatrix sample_score,
	TestType test_type)
{
	if (test_type == TestType::two_sided)
	{
		return comp_empirical_p_values_two_sided(scores, weights, sample_score);
	}
	double wei_sum = 0;
	double wei2_sum = 0;
	int sample_size = sample_score.nrow();
	if (weights.size() != sample_size)
	{
		throw std::length_error("weights has a wrong size.");
	};
	int n_scores = scores.size();
	NumericMatrix p_values(n_scores, 4);
	NumericMatrix moments(n_scores, 4);
	for (int j = 0; j < n_scores; j++)
	{
		for (int k = 0; k < 2; k++)
		{
			moments(j, k) = 0;
		}
	}
	for (int i = 0; i < sample_size; i++)
	{
		wei_sum += weights[i];
		wei2_sum += weights[i] * weights[i];
		for (int j = 0; j < n_scores; j++)
		{
			for (int k = 0; k < sample_score.ncol(); k++)
			{
				// SNP changes binding affinity
				bool stat_sign = false;
				switch (test_type)
				{
				case TestType::gte:
					stat_sign = sample_score(i, k) >= scores(j);
					break;
				case TestType::lte:
					stat_sign = sample_score(i, k) <= scores(j);
					break;
				case TestType::two_sided:
					throw std::invalid_argument("Two sided p-value should not be calculated here.");
				default:
					stat_sign = (sample_score(i, k) >= scores(j) || sample_score(i, k) <= -scores(j));
					break;
				}
				if (stat_sign)
				{
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
	for (int j = 0; j < n_scores; j++)
	{
		p_values(j, 0) = moments(j, 0) / sample_score.ncol() / sample_size;
		p_values(j, 1) = moments(j, 1) / sample_score.ncol() / sample_size - p_values(j, 0) * p_values(j, 0);
		p_values(j, 2) = p_values(j, 0) / wei_sum;
		grad2 = -p_values(j, 0) * grad1 * grad1;
		var1 = p_values(j, 1);
		cov = moments(j, 1) / sample_score.ncol() / sample_size - p_values(j, 0) * wei_sum;
		p_values(j, 1) /= sample_score.ncol() * sample_size - 1;
		if (p_values(j, 0) != wei_sum)
		{
			p_values(j, 3) = grad1 * grad1 * var1 + grad2 * grad2 * var2 + 2 * grad1 * grad2 * cov;
			p_values(j, 3) /= sample_score.ncol() * sample_size - 1;
		}
		else
		{
			p_values(j, 3) = 1;
		}
	}
	return (p_values);
}

NumericMatrix comp_empirical_p_values_two_sided(
	NumericVector scores,
	NumericVector weights,
	NumericMatrix sample_score)
{
	NumericMatrix p_values_gte = comp_empirical_p_values(scores, weights, sample_score, TestType::gte);
	NumericMatrix p_values_lte = comp_empirical_p_values(scores, weights, sample_score, TestType::lte);
	NumericMatrix p_values(scores.size(), 4);
	for (int i = 0; i < scores.size(); ++i)
	{
		for (int j = 0; j < 4; j += 2)
		{
			if (p_values_gte(i, j) > p_values_lte(i, j))
			{
				p_values(i, j) = 2 * p_values_lte(i, j);
				p_values(i, j + 1) = 4 * p_values_lte(i, j + 1);
			}
			else
			{
				p_values(i, j) = 2 * p_values_gte(i, j);
				p_values(i, j + 1) = 4 * p_values_gte(i, j + 1);
			}
		}
	}
	return p_values;
}