#include "PairedImportanceSampleClass.h"
#include "helper.h"

void PairedImportanceSamplingBase::comp_norm_const()
{
    this->norm_const = 0;
    for (int i = 0; i < cond_norm_const.size(); ++i)
    {
        this->norm_const += cond_norm_const[i];
    }
}

int PairedImportanceSamplingBase::sample_start_position(double rnd)
{
    return sample_discrete(rnd, this->cond_norm_const);
}

int PairedImportanceSamplingBase::sample_start_position()
{
    RNGScope scope;
    NumericVector rv = runif(1);
    return this->sample_start_position(rv[0]);
}

void PairedImportanceSamplingBase::initialize(double score)
{
    this->comp_cond_norm_const();
    this->comp_norm_const();
    this->set_theta(score);
}

void PairedImportanceSamplingBase::validate()
{
    if (this->theta < PairedImportanceSamplingBase::THETA_MIN)
    {
        this->theta = PairedImportanceSamplingBase::THETA_MIN;
    }
    else if (this->theta > PairedImportanceSamplingBase::THETA_MAX)
    {
        this->theta = PairedImportanceSamplingBase::THETA_MAX;
    }
    assert(this->insertion_len >= 1);
    assert(this->theta >= PairedImportanceSamplingBase::THETA_MIN && 
    this->theta <= PairedImportanceSamplingBase::THETA_MAX);
    assert(this->mc_param.stat_dist.size() == N_LETTERS);
    assert(this->mc_param.trans_mat.ncol() == N_LETTERS);
    assert(this->adj_pwm.nrow() == this->mat_d.nrow());
    assert(this->adj_pwm.ncol() == PairedImportanceSamplingBase::N_LETTERS);
    assert(this->mat_d.ncol() == PairedImportanceSamplingBase::N_LETTERS);
    rowwise_l1_normalize(mat_d, 1e-10);
}
