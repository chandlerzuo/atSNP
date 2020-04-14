#ifndef __HELPER__
#define __HELPER__

#include <Rcpp.h>
using namespace Rcpp;

int sample_discrete(double, NumericVector);
void rowwise_l1_normalize(NumericMatrix &, double);
NumericMatrix comp_empirical_p_values(
    NumericVector,
    NumericVector,
    NumericMatrix
);

#endif