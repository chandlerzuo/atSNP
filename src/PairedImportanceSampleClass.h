#ifndef __PAIRED_IS_BASE__
#define __PAIRED_IS_BASE__

#include "struct.h"
using namespace Rcpp;

// Abstract class for importance sampling
class PairedImportanceSamplingBase
{
public:
    static const double THETA_MAX, THETA_MIN;
    static const int N_LETTERS;
    MarkovChainParam mc_param;
    int insertion_len;
    double theta;
    NumericMatrix adj_pwm;
    NumericMatrix mat_d;
    NumericVector cond_norm_const;
    double norm_const;

    PairedImportanceSamplingBase(
        MarkovChainParam mc_param,
        NumericMatrix adj_pwm,
        NumericMatrix mat_d,
        int insertion_len) : mc_param(mc_param),
                             insertion_len(insertion_len),
                             theta(0),
                             adj_pwm(adj_pwm),
                             mat_d(mat_d)
    {
        this->validate();
    };
    PairedImportanceSamplingBase(
        MarkovChainParam mc_param,
        NumericMatrix adj_pwm,
        NumericMatrix mat_d,
        int insertion_len,
        double theta) : mc_param(mc_param),
                        insertion_len(insertion_len),
                        theta(theta),
                        adj_pwm(adj_pwm),
                        mat_d(mat_d)
    {
        this->validate();
    };
    int sample_start_position(double);
    int sample_start_position();
    void initialize(double);
    void comp_norm_const();
    void comp_cond_norm_const();
    virtual SampleSequence gen_importance_sample()=0;
    virtual AdjWeights gen_importance_sample_weights(IntegerVector)=0;
    virtual ScorePair comp_score_pair(NumericMatrix, IntegerVector, LoglikType)=0;

protected:
    double _comp_norm_const(NumericVector);
    
private:
    void validate();
    virtual NumericVector _comp_cond_norm_const(double)=0;
    virtual void set_theta(double)=0;
    virtual double comp_expected_score_diff(double)=0;
};

#endif