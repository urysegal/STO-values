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

static double
beta_only_f (const gsl_vector *v, void *params)
{
    auto args = static_cast<Estimator *>(params);
    return args->average_error(v);
}


void Calculated_C_Estimator::minimize(nlohmann::json &output_json)
{
    iter = 0;
    int status;

    const gsl_multimin_fdfminimizer_type *T;

    gsl_multimin_function_fdf my_func;

    my_func.n = N;
    my_func.f = beta_only_f;
    my_func.df = beta_only_df;
    my_func.fdf = beta_only_fdf;
    my_func.params = this;

    setup_initial_guess();

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    s = gsl_multimin_fdfminimizer_alloc (T, N);



    gsl_vector_view sx = gsl_vector_subvector(x, 0, N);
    gsl_multimin_fdfminimizer_set (s, &my_func, &sx.vector, step_size, tolerance);

    fprintf (stderr, "\nDirectly Calculated C:\n");

    do
    {
        update_C(s->x);

        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);
        this->estimate_error = s->f;

        if (status)
            break;

        status = gsl_multimin_test_gradient (s->gradient, stop_gradient);

        if (status == GSL_SUCCESS)
            fprintf (stderr, "Minimum found at:\n");

        fprintf (stderr,"%5lu beta=%.5f C=%.5f f=%10.5f\n", iter,
                 gsl_vector_get (s->x, 0),
                 gsl_vector_get (x, 1),
                 s->f);

    }
    while (status == GSL_CONTINUE && iter < args.get_max_iterations());

    output_results(output_json, x, s->x);
}

}