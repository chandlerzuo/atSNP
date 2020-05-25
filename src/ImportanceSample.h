#ifndef __IS__
#define __IS__

#include "MotifScore.h"

NumericMatrix gen_utility_matrix(NumericMatrix, NumericMatrix, int, double);
NumericVector gen_prob_start_pos(NumericMatrix, int, NumericVector);
RcppExport SEXP test_gen_utility_matrix(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP test_gen_prob_start_pos(SEXP, SEXP, SEXP);

#endif