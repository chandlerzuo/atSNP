#ifndef __IS_INDEL__
#define __IS_INDEL__

#include "PairedImportanceSampleClass.h"
using namespace Rcpp;

class ImportanceSampleIndel : public PairedImportanceSamplingBase
{
public:
    ImportanceSampleIndel(
        MarkovChainParam mc_param,
        NumericMatrix adj_pwm,
        NumericMatrix mat_d,
        int insertion_len) : PairedImportanceSamplingBase(mc_param,
                                                          adj_pwm,
                                                          mat_d,
                                                          insertion_len){};
    ImportanceSampleIndel(
        MarkovChainParam mc_param,
        NumericMatrix adj_pwm,
        NumericMatrix mat_d,
        int insertion_len,
        double theta) : PairedImportanceSamplingBase(mc_param,
                                                     adj_pwm,
                                                     mat_d,
                                                     insertion_len,
                                                     theta){};
    SampleSequence gen_importance_sample();
    AdjWeights gen_importance_sample_weights(IntegerVector);
    ScorePair comp_score_pair(NumericMatrix, IntegerVector, LoglikType);

private:
    NumericVector _comp_cond_norm_const(double);
    void set_theta(double);
    double comp_expected_score_diff(double);
};

RcppExport SEXP p_value_change_indel(
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP
); 

RcppExport SEXP comp_indel_motif_scores(SEXP, SEXP, SEXP);

RcppExport SEXP test_importance_sample_indel(
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP,
    SEXP);

#endif