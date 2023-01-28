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

#define GET_beta(i) gsl_vector_get(v, (i))

using namespace std;

namespace stovalues {


real_t pi = std::numbers::pi_v<long double>;
real_t sqrt_pi = sqrtq(pi);
real_t half_sqrt_pi = sqrtq(pi)/2.0;
real_t quarter_sqrt_pi = sqrtq(pi)/4.0;


real_t Estimator::errfunc(real_t sqrt_beta)
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
        sum += Cj / (sqrtq(beta_i + beta_j));
    }
    sum *= sqrt_pi;

    real_t sqrt_beta_i = sqrtq(beta_i);
    real_t term = (1.0/sqrtq(beta_i))*expq(1/(4.0*beta_i)) * errfunc(sqrt_beta_i);
    sum -= sqrt_pi*term;

    return sum;
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
        sum += Ci*Cj*1.0/sqrtq(powq(beta_i+beta_j,3.0)) ;
    }
    sum *= -half_sqrt_pi;

    auto err = errfunc(sqrtq(beta_i));
    auto exp = expq(1/(4*beta_i));

    sum += half_sqrt_pi * Ci * beta_i_pow_neg_3_2* exp * err ;
    sum += quarter_sqrt_pi * Ci * beta_i_pow_neg_5_2 * exp * err ;
    sum -= Ci/(2*beta_i*beta_i);

    return sum;
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

static double
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



Estimator::~Estimator()
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


void Estimator::setup_initial_guess()
{
    x = gsl_vector_alloc (N*2);
    if ( N ==1 ) {
        gsl_vector_set(x, 0, 0.025);
        gsl_vector_set(x, 1, 1);
    } else {
        try {
            std::ifstream ifs(args.get_guess_file());
            nlohmann::json jf = nlohmann::json::parse(ifs);
            if ( jf["N"] >= N-1 ) {
                logger()->critical("Guess file {} does not contain enough terms: has {}. {} needed.",
                                   args.get_guess_file(), jf["N"], N-1);
                exit(1);
            }
#if 0
            double dat[] = {
            0.3189040000E+01     ,  0.3627800244E+00,  0.9348630000E+01    ,   0.4684630315E+00
            };
            for ( auto i = 0U ; i < 4 ; ++i ) {
               gsl_vector_set(x, i, dat[i]);
//               gsl_vector_set(x, i+N, (i+1)*2);
            }
#else
            string best_method = jf["best_method"];
            nlohmann::json method_output = jf[best_method];
            nlohmann::json terms_output = method_output["terms"];
            unsigned int i = 0;
            double sum_beta = 0;
            double last_beta = 0;
            double last_C = 0;
            for ( auto &it: terms_output ) {
                if ( i == N-1) {
                    break;
                }
                last_beta = double(it["beta"]);
                sum_beta += last_beta;
                gsl_vector_set(x, i, last_beta);
                last_C = it["C"];
                gsl_vector_set(x, i+N, last_C);
                ++i;
            }
            gsl_vector_set(x, N-1, sum_beta*1.1);
            //this->update_C(x);
            gsl_vector_set(x, (2*N)-1, last_C/double (N));
#endif
        } catch (std::exception &e) {
            logger()->critical("Reading Guess: FATAL: {}", std::string(e.what()));
            exit(1);
        }

    }
}

void Estimator::output_results(nlohmann::json &output_json, const gsl_vector *C_vector, const gsl_vector *beta_vector)
{
    double test_x = 0.1;
    double real_res = exp(-test_x);
    double estimate = 0;
    for ( auto i = 0U ; i < N ; i++ ) {
        auto C = gsl_vector_get(C_vector, N+i) ;
        auto beta = gsl_vector_get(beta_vector, i);
        estimate += C* exp(-beta*(test_x*test_x));
        result_term new_term= { gsl_vector_get(C_vector, N+i), gsl_vector_get(beta_vector, i)} ;
        this->terms.emplace_back(new_term);
    }
    fprintf(stderr, "At X=%f, %f %f (err = %f)\n", test_x, real_res, estimate, abs(real_res-estimate));
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
    std::cerr << result << std::endl;
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

    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tolerance);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        this->estimate_error = s->f;

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, stop_gradient);

        if (status == GSL_SUCCESS)
            fprintf (stderr,"Minimum found at:\n");

        if ( iter % 100 == 1 ) {
            fprintf(stderr, "%u %5lu f=%10.15f\n", N, iter,
                    s->f);
        }
    }
    while (status == GSL_CONTINUE && iter < args.get_max_iterations());
    output_results(output_json, s->x, s->x);


}

} // namespace


