#include <iostream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>

#include <numbers>
#include <quadmath.h>

#include "s2g.h"
#include "logger.h"


using namespace std;

namespace stovalues {


static real_t pi = std::numbers::pi_v<long double>;
static real_t sqrt_pi = sqrtq(pi);
static real_t half_sqrt_pi = sqrtq(pi) / 2.0;
static real_t quarter_sqrt_pi = sqrtq(pi) / 4.0;

Three_D_Estimator::~Three_D_Estimator()
{

}

real_t
Three_D_Estimator::get_Bi(const gsl_vector *v, int i) const
{
    real_t beta_i = GET_beta(v, i);

    real_t term = ( (2*beta_i) +1 ) / sqrtq(powq(beta_i, 5)) ;
    term *= expq(1.0/(4.0*beta_i));
    real_t err= errfunc(sqrtq(beta_i));
    return term*err;
}



void Three_D_Estimator::adjust_C(gsl_vector *v)
{
    auto C_adjuster = this->calculate_adjuster(v);
    for ( auto j= 0U ; j < N ; j++ ) {
        gsl_vector_set(v, N + j, GET_C(v, j) * C_adjuster);
    }
}

// STO-2G : C of 1 = 0.208388432054940 , beta = .151623
//               2 = .481774983399508 , beta = .851819

double Three_D_Estimator::average_error (const gsl_vector *original_v)
{
    real_t double_sum=0, square_sum=0, B_sum = 0;

    auto v = gsl_vector_alloc(N*2);
    for ( auto k = 0U ; k < 2*N ; ++k ) {
        gsl_vector_set(v, k, gsl_vector_get(original_v, k));
    }
    //adjust_C(v);

    for ( auto i=0U ; i< N ; ++i ) {

        auto Ci = GET_C(v, i);
        auto beta_i = GET_beta(v, i);
        auto Bi = get_Bi(v, i);

        B_sum += Ci*Bi;
        square_sum += Ci / (2.0*beta_i*beta_i);

        for ( auto j = 0U ; j < N ; ++j ) {
            auto Cj = GET_C(v, j);
            auto beta_j = GET_beta(v, j);
            double_sum += Ci*Cj/(sqrtq(powq(beta_i+beta_j,3)));
        }

    }
    //this->C_adjuster = 1.0/(sqrtq((double_sum*sqrt_pi)));
    //printf("Adjuster: %10.15f, N=%10.15f\n", double(C_adjuster), double(double_sum*sqrt_pi));

    gsl_vector_free(v);
    return 1.0/4.0 + quarter_sqrt_pi*double_sum + square_sum - quarter_sqrt_pi*B_sum;

}

double Three_D_Estimator::diff_by_Ci(const gsl_vector *original_v, size_t i)
{
    auto v = gsl_vector_alloc(N*2);
    for ( auto k = 0U ; k < 2*N ; ++k ) {
        gsl_vector_set(v, k, gsl_vector_get(original_v, k));
    }
   // adjust_C(v);

    real_t beta_i = GET_beta(v, i);
    auto Bi= get_Bi(v, i);
    real_t sum = 0;
    for ( auto j = 0U ; j < N ; ++j )
    {
        real_t beta_j = GET_beta(v, j);
        sum += GET_C(v, j) / ( sqrtq(powq(beta_i+beta_j,3)));
    }
    gsl_vector_free(v);

    return half_sqrt_pi * sum + (1.0/(2*beta_i*beta_i)) - quarter_sqrt_pi*Bi;
}


double Three_D_Estimator::diff_by_bi(const gsl_vector *original_v, size_t i)
{

    auto v = gsl_vector_alloc(N*2);
    for ( auto k = 0U ; k < 2*N ; ++k ) {
        gsl_vector_set(v, k, gsl_vector_get(original_v, k));
    }
   // adjust_C(v);

    auto Ci = GET_C(v, i);
    auto beta_i = GET_beta(v, i);
    auto Bi = get_Bi(v, i);



    real_t sum = 0;
    for ( auto j = 0U ; j < N ; ++j )
    {
        real_t beta_j = GET_beta(v, j);
        real_t Cj = GET_C(v, j);
        sum +=  Ci*Cj/sqrtq(powq(beta_i+beta_j,5));
    }

    real_t frac = 12.0*beta_i*beta_i + 12*beta_i + 1;
    frac /= 16.0*beta_i*beta_i*(2.0*beta_i + 1);

    real_t deriv_term = sqrt_pi * frac * Bi ;
    deriv_term -= (2.0*beta_i + 1) / (8.0*powq(beta_i, 4)) ;
    deriv_term *= Ci;
    gsl_vector_free(v);

    return -3*quarter_sqrt_pi*sum - Ci/powq(beta_i, 3) + deriv_term;
}



real_t Three_D_Estimator::calculate_adjuster (const gsl_vector *v)
{
    real_t double_sum=0;
    real_t C_adjuster = 0;

    for ( auto i=0U ; i< N ; ++i ) {

        auto Ci = GET_C(v, i);
        auto beta_i = GET_beta(v, i);

        for ( auto j = 0U ; j < N ; ++j ) {
            auto Cj = GET_C(v, j);
            auto beta_j = GET_beta(v, j);
            double_sum += Ci*Cj/(sqrtq(powq(beta_i+beta_j,3)));
        }

    }
    C_adjuster = 1.0/(sqrtq((double_sum*sqrt_pi)));
    //printf("Adjuster: %10.15f, N=%10.15f\n", double(C_adjuster), double(double_sum*sqrt_pi));

    return C_adjuster;
}


}



