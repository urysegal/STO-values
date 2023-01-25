#include "s2g.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>
#include "logger.h"
#include <numbers>
#include <quadmath.h>


namespace stovalues {

typedef  __float128 real_t;

real_t pi = std::numbers::pi_v<long double>;
real_t half_pi = std::numbers::pi_v<long double> / 2.0 ;
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

    for ( auto i = 0U ; i < v->size ; i+=2 ) {

        real_t Ci = gsl_vector_get(v, i);
        real_t beta_i = gsl_vector_get(v, i+1);

        real_t sqrt_beta_i = sqrtq(beta_i);
        sum2 += (Ci/sqrt_beta_i)*expq(1.0/(4.0*beta_i))*errfunc(sqrt_beta_i);

        for ( auto j = 0U ; j < v->size ; j+=2 ) {

            real_t Cj = gsl_vector_get(v, j);
            auto beta_j = gsl_vector_get(v, j + 1);

            sum1 += Ci * Cj / (sqrtq(beta_i + beta_j));
        }
    }

    double sum = 0.5 + half_pi*sum1 + pi*sum2;

    logger()->trace("f called with N={}, sum is {}",v->size/2,sum);

    return sum;
}

double Estimator::diff_by_Ci(const gsl_vector *v, size_t i)
{

    real_t beta_i = gsl_vector_get(v, i + 1);
    real_t sum = 0;
    for (auto j = 0U; j < v->size; j += 2) {

        real_t Cj = gsl_vector_get(v, j);
        real_t beta_j = gsl_vector_get(v, j + 1);
        sum += Cj/(sqrtq(beta_i+beta_j));

        real_t sqrt_beta_j = sqrtq(beta_j);
        real_t term = expq(1/(4.0*beta_j)) * errfunc(sqrt_beta_j);
        sum -= sqrt_pi*term;
    }

    return sqrt_pi*sum;
}

double Estimator::diff_by_bi(const gsl_vector *v, size_t i)
{
    real_t beta_i = gsl_vector_get(v, i + 1);
    real_t beta_i_sqrt = sqrtq(beta_i);
    real_t beta_i_pow_neg_3_2 = 1.0/powq(beta_i_sqrt, 3.0);
    real_t beta_i_pow_neg_5_2 = 1.0/powq(beta_i_sqrt, 5.0);
    real_t Ci = gsl_vector_get(v, i);

    real_t sum = 0;
    for (auto j = 0U; j < v->size; j += 2) {

        real_t Cj = gsl_vector_get(v, j);
        real_t beta_j = gsl_vector_get(v, j + 1);

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
    auto N = v->size/2;
    logger()->trace("df called with N={}", N);

    for ( auto i = 0U ; i < v->size ; i+=2 ) {
        auto dF_dCi = diff_by_Ci(v, i);
        gsl_vector_set(df, i, dF_dCi);
        auto dF_dbi = diff_by_bi(v, i);
        gsl_vector_set(df, i+1, dF_dbi);
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

void Estimator::work(nlohmann::json &output_json)
{
    size_t iter = 0;
    int status;
    auto N = args.get_number_of_terms();
    int n = N*2;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;

    gsl_vector *x;
    gsl_multimin_function_fdf my_func;

    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = this;

    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, 1.0);
    gsl_vector_set (x, 1, 0.27);

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, n);

    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, 1e-3);

        if (status == GSL_SUCCESS)
            printf ("Minimum found at:\n");

        printf ("%5lu %.5f %.5f %10.5f\n", iter,
                gsl_vector_get (s->x, 0),
                gsl_vector_get (s->x, 1),
                s->f);

    }
    while (status == GSL_CONTINUE && iter < args.get_max_iterations());

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

}

} // namespace


