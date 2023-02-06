#include <iostream>
#include <fstream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>

#include <numbers>
#include <quadmath.h>

#include "s2g.h"
#include "logger.h"


using namespace std;

namespace stovalues {


real_t pi = std::numbers::pi_v<long double>;
real_t sqrt_pi = sqrtq(pi);
real_t half_sqrt_pi = sqrtq(pi)/2.0;
real_t quarter_sqrt_pi = sqrtq(pi)/4.0;


real_t Guess_Estimator::errfunc(real_t sqrt_beta)
{
    return 1.0 - gsl_sf_erf(1.0/(2.0*sqrt_beta));
}


double Guess_Estimator::average_error(const gsl_vector *v)
{
    real_t sum1 = 0;
    real_t sum2 = 0;

    for ( auto i = 0U ; i < N ; ++i ) {

        real_t Ci = GET_C(v, i);
        real_t beta_i = GET_beta(v, i);

        real_t sqrt_beta_i = sqrtq(beta_i);
        sum2 += (Ci/sqrt_beta_i)*expq(1.0/(4.0*beta_i))*errfunc(sqrt_beta_i);

        for ( auto j = 0U ; j < N ; ++j ) {

            real_t Cj = GET_C(v, j);
            auto beta_j = GET_beta(v, j);

            sum1 += Ci * Cj / (sqrtq(beta_i + beta_j));
        }
    }

    double sum = 0.5 + half_sqrt_pi*sum1 - sqrt_pi*sum2;

    //logger()->trace("f called with N={}, sum is {}",v->size/2,sum);

    return sum;
}

double Guess_Estimator::diff_by_Ci(const gsl_vector *v, size_t i)
{

    real_t beta_i = GET_beta(v, i);
    real_t sum = 0;
    for (auto j = 0U; j < N; ++j) {

        real_t Cj = GET_C(v, j);
        real_t beta_j = GET_beta(v, j);
        sum += Cj / (sqrtq(beta_i + beta_j));
    }
    sum *= sqrt_pi;

    real_t sqrt_beta_i = sqrtq(beta_i);
    real_t term = (1.0/sqrtq(beta_i))*expq(1/(4.0*beta_i)) * errfunc(sqrt_beta_i);
    sum -= sqrt_pi*term;

    return sum;
}

double Guess_Estimator::diff_by_bi(const gsl_vector *v, size_t i)
{
    real_t beta_i = GET_beta(v, i);
    real_t beta_i_sqrt = sqrtq(beta_i);
    real_t beta_i_pow_neg_3_2 = 1.0/powq(beta_i_sqrt, 3.0);
    real_t beta_i_pow_neg_5_2 = 1.0/powq(beta_i_sqrt, 5.0);
    real_t Ci = GET_C(v, i);

    real_t sum = 0;
    for (auto j = 0U; j < N; ++j) {
        real_t Cj = GET_C(v, j);
        real_t beta_j = GET_beta(v, j);
        sum += Ci*Cj*1.0/sqrtq(powq(beta_i+beta_j,3.0)) ;
    }
    sum *= -half_sqrt_pi;

    auto err = errfunc(beta_i_sqrt);
    auto exp = expq(1/(4*beta_i));

    sum += half_sqrt_pi * Ci * beta_i_pow_neg_3_2* exp * err ;
    sum += quarter_sqrt_pi * Ci * beta_i_pow_neg_5_2 * exp * err ;
    sum -= Ci/(2*beta_i*beta_i);

    return sum;
}


void Guess_Estimator::average_error_df (const gsl_vector *v, gsl_vector *df)
{
    for ( auto i = 0U ; i < N ; ++i ) {
        auto dF_dCi = diff_by_Ci(v, i);
        gsl_vector_set(df, N+i, dF_dCi);
        auto dF_dbi = diff_by_bi(v, i);
        gsl_vector_set(df, i, dF_dbi);
       // logger()->trace("dF/dC_{} = {} , dF/db_{} = {}", i, dF_dCi, i, dF_dbi);
    }
}



/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    auto args = static_cast<Guess_Estimator *>(params);

    *f = args->average_error(x);
    args->average_error_df(x, df);
}

static double
my_f (const gsl_vector *v, void *params)
{
    auto args = static_cast<Guess_Estimator *>(params);
    return args->average_error(v);
}

void
my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
    auto args = static_cast<Guess_Estimator *>(params);
    args->average_error_df(v,df);
}



Guess_Estimator::~Guess_Estimator()
{
    gsl_vector_free (x);
    x=nullptr;
    gsl_multimin_fdfminimizer_free (s);
    s = nullptr;
}




void Guess_Estimator::setup_initial_guess()
{
    x = gsl_vector_alloc (N*2);
    if ( N ==1 ) {
        gsl_vector_set(x, 0, 0.025);
        gsl_vector_set(x, 1, 1);
    } else {
        assert (args.get_max_guesses()) ;
        double beta = args.get_initial_beta();
        if ( args.get_max_guesses() > 1 ) {
            beta += drand48() * 0.0001;
        }
        for (auto i = 0U; i < N; ++i) {
            gsl_vector_set(x, i, beta);
            if ( args.get_max_guesses() > 1 ) {
                beta *= (2 - drand48() * 0.99);
            } else {
                beta *= args.get_beta_factor();
            }
            if ( args.get_max_guesses() > 1 ) {
                gsl_vector_set(x, i + N, args.get_initial_C() + drand48() * 0.1);
            } else {
                gsl_vector_set(x, i + N, args.get_initial_C());
            }
        }
    }
}

void Guess_Estimator::output_results(nlohmann::json &output_json, const gsl_vector *C_vector, const gsl_vector *beta_vector)
{

    for ( auto i = 0U ; i < N ; i++ ) {
        auto C = gsl_vector_get(C_vector, N+i) ;
        auto beta = gsl_vector_get(beta_vector, i);
        result_term new_term= { C, beta} ;
        this->terms.emplace_back(new_term);
    }

    std::vector<double> matlab_C, matlab_beta;
    output_set["error"] = this->estimate_error;
    output_set["iterations"] = iter;
    std::vector<nlohmann::json> json_terms;
    for ( auto const &it : terms ) {
        nlohmann::json term;
        term["C"] = it.C;
        term["beta"] = it.beta;
        matlab_C.emplace_back(it.C);
        matlab_beta.emplace_back(it.beta);
        json_terms.emplace_back(term);
    }
    output_set["terms"] = json_terms;
    output_set["matlab_C"] = matlab_C;
    output_set["matlab_beta"] = matlab_beta;

}





void Guess_Estimator::minimize(nlohmann::json &output_json)
{
    iter = 0;
    int status;
    int n = N*2;

    srand48(::time(nullptr));

    const gsl_multimin_fdfminimizer_type *T;

    gsl_multimin_function_fdf my_func;

    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = this;

    setup_initial_guess();

    auto last_good_result = gsl_vector_alloc (N*2);

    T = gsl_multimin_fdfminimizer_conjugate_pr;
    s = gsl_multimin_fdfminimizer_alloc (T, n);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tolerance);


    this->estimate_error = 1e10;
    do
    {
        iter++;
        gsl_vector_memcpy(last_good_result, x);

        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, stop_gradient);


        if ( isnan(s->f)) {
            break;
        } else {
            this->estimate_error = s->f;
        }

    }
    while ( status == GSL_CONTINUE && iter < args.get_max_iterations());
    if (isnan(s->f)) {
        gsl_vector_memcpy(s->x, last_good_result);
    }
    output_results(output_json, s->x, s->x);
    gsl_vector_free(last_good_result);

}

} // namespace


