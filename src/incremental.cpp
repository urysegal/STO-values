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

Incremental_Estimator::~Incremental_Estimator()
{
    gsl_vector_free (adjusted_x);
    adjusted_x=nullptr;
}


void Incremental_Estimator::setup_initial_guess()
{
    x = gsl_vector_alloc (N*2);
    adjusted_x = gsl_vector_alloc (N*2);

    if ( N ==1 ) {
        assert("Cannot increament result with zero previous results. Use a guess");
    } else {
        nlohmann::json terms_output;

        try {
            std::ifstream ifs(args.get_guess_file());
            nlohmann::json jf = nlohmann::json::parse(ifs);
            if (jf["N"] >= N - 1) {
                logger()->critical("Guess file {} does not contain enough terms: has {}. {} needed.",
                                   args.get_guess_file(), jf["N"], N - 1);
                exit(1);
            }

            string best_method = jf["best_method"];
            nlohmann::json method_output = jf[best_method];
            terms_output = method_output["terms"];
        } catch (std::exception &e) {
            logger()->critical("Reading Guess: FATAL: {}", std::string(e.what()));
            exit(1);
        }

        unsigned int i = 0;
        double last_beta = 0;
        double last_C = 0;
        for (auto &it: terms_output) {
            if (i == N - 1) {
                break;
            }
            last_beta = double(it["beta"]);
            gsl_vector_set(x, i, last_beta);
            last_C = it["C"];
            gsl_vector_set(x, i + N, last_C);
            ++i;
        }
    }
}


static void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    auto args = static_cast<Incremental_Estimator *>(params);

    *f = args->average_error(x);
    args->average_error_df(x, df);
}

static double
my_f (const gsl_vector *v, void *params)
{
    auto args = static_cast<Incremental_Estimator *>(params);
    return args->average_error(v);
}

static void
my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
    auto args = static_cast<Incremental_Estimator *>(params);
    args->average_error_df(v,df);
}

void Incremental_Estimator::calculate_adjusted_x(const gsl_vector *v)
{
    // adjust beta - old betas
    for ( auto i = 0U ; i < N-1;  ++i ) {
        gsl_vector_set(adjusted_x, i, GET_beta(x,i));
    }
    // and one new beta
    gsl_vector_set(adjusted_x, N-1, gsl_vector_get(v, 0));

    // adjust C - old C, adjusted by factor
    auto factor = gsl_vector_get(v, 2);
    for ( auto i = 0U ; i < N-1;  ++i ) {
        gsl_vector_set(adjusted_x, N+i, factor * GET_C(x,i));
    }
    // and one new C
    gsl_vector_set(adjusted_x, (N*2)-1, gsl_vector_get(v, 1));

}

double Incremental_Estimator::average_error(const gsl_vector *v)
{
    calculate_adjusted_x(v);

    return Guess_Estimator::average_error(adjusted_x);
}

void Incremental_Estimator::average_error_df(const gsl_vector *v, gsl_vector *df)
{
    calculate_adjusted_x(v);

    auto dF_dCi = diff_by_Ci(adjusted_x, N-1);
    auto dF_dbi = diff_by_bi(adjusted_x, N-1);

    gsl_vector_set(df, 0, dF_dbi);
    gsl_vector_set(df, 1, dF_dCi);

    double dF_factor = 0;
    for ( auto i = 0U ; i < N-1 ; ++i ) {
        auto dCi = diff_by_Ci(adjusted_x, i);
        dF_factor += dCi;
    }
    gsl_vector_set(df, 2, dF_factor);

}

void Incremental_Estimator::minimize(nlohmann::json &output_json)
{
    iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;

    gsl_multimin_function_fdf my_func;

    my_func.n = 3; // One value controls the coefficient on the previous result, one new C one new beta
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = this;

    setup_initial_guess();

    T = gsl_multimin_fdfminimizer_conjugate_pr;
    s = gsl_multimin_fdfminimizer_alloc (T, 3);

    auto new_params = gsl_vector_alloc(3);
    gsl_vector_set(new_params,0, GET_beta(x, N-2)*1.01);
    gsl_vector_set(new_params,1,1/N);
    gsl_vector_set(new_params,2,0.8); // Guess coefficient for previous result


    gsl_multimin_fdfminimizer_set (s, &my_func, new_params, step_size, tolerance);

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, stop_gradient);


        //if (status == GSL_SUCCESS)
        //  fprintf (stderr,"Minimum found at:\n");

        // if ( iter % 100 == 1 ) {
        fprintf(stderr, "%u %5lu f=%10.15f beta=%10.15f C=%10.15f factor=%10.15f \n", N, iter,
                 s->f, gsl_vector_get(s->x,0),gsl_vector_get(s->x,1), gsl_vector_get(s->x,2) );
        //  }

        if ( isnan(s->f)) {
            this->estimate_error = 9999999999;
        } else {
            this->estimate_error = s->f;
        }

    }
    while (not isnan(s->f) and  status == GSL_CONTINUE && iter < args.get_max_iterations());

    calculate_adjusted_x(s->x);
    output_results(output_json, adjusted_x, adjusted_x);

    gsl_vector_free(new_params);


}
}