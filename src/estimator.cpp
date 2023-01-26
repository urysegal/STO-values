#include <iostream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_linalg.h>

#include <numbers>
#include <quadmath.h>

#include "s2g.h"
#include "logger.h"

#define GET_beta(i) gsl_vector_get(v, (i))

namespace stovalues {

typedef  __float128 real_t;

real_t pi = std::numbers::pi_v<long double>;
real_t sqrt_pi = sqrtq(pi);
real_t half_sqrt_pi = sqrtq(pi)/2.0;
real_t quarter_sqrt_pi = sqrtq(pi)/4.0;


real_t errfunc(real_t sqrt_beta)
{
    return 1.0 - gsl_sf_erf(1.0/(2.0*sqrt_beta));
}


double Estimator::average_error(const gsl_vector *v)
{
    real_t sum1 = 0;
    real_t sum2 = 0;

    for ( auto i = 0U ; i < N ; ++i ) {

        real_t Ci = GET_C(v, i);
        real_t beta_i = GET_beta(i);

        real_t sqrt_beta_i = sqrtq(beta_i);
        sum2 += (Ci/sqrt_beta_i)*expq(1.0/(4.0*beta_i))*errfunc(sqrt_beta_i);

        for ( auto j = 0U ; j < N ; ++j ) {

            real_t Cj = GET_C(v, j);
            auto beta_j = GET_beta(j);

            sum1 += Ci * Cj / (sqrtq(beta_i + beta_j));
        }
    }

    double sum = 0.5 + half_sqrt_pi*sum1 - sqrt_pi*sum2;

    logger()->trace("f called with N={}, sum is {}",v->size/2,sum);

    return sum;
}

double Estimator::diff_by_Ci(const gsl_vector *v, size_t i)
{

    real_t beta_i = GET_beta(i);
    real_t sum = 0;
    for (auto j = 0U; j < N; ++j) {

        real_t Cj = GET_C(v, j);
        real_t beta_j = GET_beta(j);
        sum += Cj/(sqrtq(beta_i+beta_j));

        real_t sqrt_beta_j = sqrtq(beta_j);
        real_t term = (1.0/sqrtq(beta_j))*expq(1/(4.0*beta_j)) * errfunc(sqrt_beta_j);
        sum -= sqrt_pi*term;
    }

    return sqrt_pi*sum;
}

double Estimator::diff_by_bi(const gsl_vector *v, size_t i)
{
    real_t beta_i = GET_beta(i);
    real_t beta_i_sqrt = sqrtq(beta_i);
    real_t beta_i_pow_neg_3_2 = 1.0/powq(beta_i_sqrt, 3.0);
    real_t beta_i_pow_neg_5_2 = 1.0/powq(beta_i_sqrt, 5.0);
    real_t Ci = GET_C(v, i);

    real_t sum = 0;
    for (auto j = 0U; j < N; ++j) {

        real_t Cj = GET_C(v, j);
        real_t beta_j = GET_beta(j);

        auto err = errfunc(sqrtq(beta_j));
        auto exp = expq(1/(4*beta_j));
        sum += Ci*Cj*1.0/sqrtq(powq(beta_i+beta_j,3.0)) ;
        sum += half_sqrt_pi * Ci * beta_i_pow_neg_3_2* exp * err ;
        sum += quarter_sqrt_pi * Ci * beta_i_pow_neg_5_2 * exp * err ;
        sum -= Ci/(2*beta_j*beta_j);
    }
    return -half_sqrt_pi*sum;
}


void Estimator::average_error_df (const gsl_vector *v, gsl_vector *df)
{
    for ( auto i = 0U ; i < N ; ++i ) {
        auto dF_dCi = diff_by_Ci(v, i);
        gsl_vector_set(df, N+i, dF_dCi);
        auto dF_dbi = diff_by_bi(v, i);
        gsl_vector_set(df, i, dF_dbi);
        logger()->trace("dF/dC_{} = {} , dF/db_{} = {}", i, dF_dCi, i, dF_dbi);
    }
}



/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    auto args = static_cast<Estimator *>(params);

    *f = args->average_error(x);
    args->average_error_df(x, df);
}

double
my_f (const gsl_vector *v, void *params)
{
    auto args = static_cast<Estimator *>(params);
    return args->average_error(v);
}

void
my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
    auto args = static_cast<Estimator *>(params);
    args->average_error_df(v,df);
}

void Estimator::minimize(nlohmann::json &output_json)
{
    iter = 0;
    int status;
    int n = N*2;

    const gsl_multimin_fdfminimizer_type *T;

    gsl_multimin_function_fdf my_func;

    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = this;

    setup_initial_guess();

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, n);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 0.1);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        this->estimate_error = s->f;

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS)
            fprintf (stderr,"Minimum found at:\n");

        fprintf (stderr,"%5lu beta=%.5f C=%.5f f=%10.5f\n", iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->f);

    }
    while (status == GSL_CONTINUE && iter < args.get_max_iterations());
    output_results(output_json, s->x, s->x);


}

Estimator::~Estimator() noexcept
{
    gsl_vector_free (x);
    x=nullptr;
    gsl_multimin_fdfminimizer_free (s);
    s = nullptr;
}


void Calculated_C_Estimator::average_error_df_beta_only(const gsl_vector *v, gsl_vector *df)
{
    for ( auto i = 0U ; i < N ; ++i ) {
        auto dF_dbi = diff_by_bi(v, i);
        gsl_vector_set(df, i, dF_dbi);
        logger()->trace("dF/db_{} = {}", i, dF_dbi);
    }
}

/* Compute both f and df together. */
void
beta_only_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    auto args = static_cast<Calculated_C_Estimator *>(params);

    *f = args->average_error(x);
    args->average_error_df_beta_only(x, df);
}

void
beta_only_df (const gsl_vector *v, void *params, gsl_vector *df)
{
    auto args = static_cast<Calculated_C_Estimator *>(params);
    args->average_error_df_beta_only(v,df);
}

void Estimator::setup_initial_guess()
{
    x = gsl_vector_alloc (N*2);
    if ( N ==1 ) {
        gsl_vector_set(x, 0, 0.27);
        gsl_vector_set(x, 1, 1);
    } else {

    }
}

void Estimator::output_results(nlohmann::json &output_json, const gsl_vector *C_vector, const gsl_vector *beta_vector)
{
    for ( auto i = 0U ; i < N ; i++ ) {
        result_term new_term= { gsl_vector_get(C_vector, N+i), gsl_vector_get(beta_vector, i)} ;
        this->terms.emplace_back(new_term);
    }

    nlohmann::json result;
    result["error"] = this->estimate_error;
    result["iterations"] = iter;
    std::vector<nlohmann::json> json_terms;
    for ( auto const &it : terms ) {
        nlohmann::json term;
        term["C"] = it.C;
        term["beta"] = it.beta;
        json_terms.emplace_back(term);
    }
    result["terms"] = json_terms;
    output_json[get_method_name()] = result;
    std::cout << result << std::endl;
}

void Calculated_C_Estimator::minimize(nlohmann::json &output_json)
{
    iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;

    gsl_multimin_function_fdf my_func;

    my_func.n = N;
    my_func.f = my_f;
    my_func.df = beta_only_df;
    my_func.fdf = beta_only_fdf;
    my_func.params = this;

    setup_initial_guess();

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, N);

    gsl_vector_view sx = gsl_vector_subvector(x, 0, N);
    gsl_multimin_fdfminimizer_set (s, &my_func, &sx.vector, 0.01, 0.1);

    printf ("\nDirectly Calculated C:\n");

    do
    {
        update_C(s->x);

        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);
        this->estimate_error = s->f;

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        fprintf (stderr,"%5lu beta=%.5f C=%.5f f=%10.5f\n", iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (x, 1),
                s->f);

    }
    while (status == GSL_CONTINUE && iter < args.get_max_iterations());

    output_results(output_json, x, s->x);

}

void Estimator::update_C(gsl_vector *betas)
{
    gsl_matrix *Ski = gsl_matrix_alloc(N, N);
    gsl_vector *S0i = gsl_vector_alloc(N);
    gsl_vector *C = gsl_vector_alloc (N);

    gsl_permutation * p = gsl_permutation_alloc (N);

    for ( auto k = 0U ; k < N ; ++k ) {
        for ( auto i = 0U ; i < N ; ++i ) {
            real_t beta_k = gsl_vector_get(betas, k);
            real_t beta_i = gsl_vector_get(betas, i);
            auto term = 1.0/sqrtq(beta_i+beta_k);
            gsl_matrix_set(Ski,k,i,term);
        }
    }

    for ( auto i = 0U ; i < N ; ++i ) {
        real_t beta_i = gsl_vector_get(betas, i);
        real_t sqrt_beta_i = sqrtq(beta_i) ;
        auto erf = errfunc(sqrt_beta_i) ;
        auto term = (1/ sqrt_beta_i)*expq(1.0/4.0*beta_i)* erf;
        gsl_vector_set(S0i, i, term);
    }

    int decomp_sign = 0;
    int status = gsl_linalg_LU_decomp (Ski, p, &decomp_sign);
    if ( status != GSL_SUCCESS ) {
        throw std::runtime_error("LU Decomposition failed");
    }

    status = gsl_linalg_LU_solve (Ski, p, S0i, C);
    if ( status != GSL_SUCCESS ) {
        throw std::runtime_error("LU matrix solver failed");
    }

//    printf ("New C = \n");
//    gsl_vector_fprintf (stdout, C, "%g");
    for ( auto i = 0U ; i < N ; ++i ) {
        gsl_vector_set(x,i+N, gsl_vector_get(C,i) );
    }

    gsl_permutation_free (p);
    gsl_matrix_free(Ski);
    gsl_vector_free(S0i);
    gsl_vector_free(C);

}

} // namespace


