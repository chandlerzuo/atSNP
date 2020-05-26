#ifndef __STRUCT__
#define __STRUCT__

#include <Rcpp.h>
using namespace Rcpp;

enum LoglikType {max, mean, median};
enum TestType {gte,lte,absolute,two_sided};

struct SequenceScores {
  int best_match_pos;
  float max_log_lik;
  float mean_log_lik;
  float median_log_lik;
};

// Markov Chain parameter
struct MarkovChainParam
{
    NumericVector stat_dist; // stationary distribution
    NumericMatrix trans_mat; // transition matrix
};

// Sample Results
struct SampleSequence
{
    int start_pos;
    IntegerVector sequence;
};

// Adjustment Ratios
struct AdjWeights
{
    double joint;
    double base;
};

// Scores
struct ScorePair
{
    double long_seq_score;
    double short_seq_score;
};

#endif