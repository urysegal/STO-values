#include "s2g.h"
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include "logger.h"
#include <numbers>


namespace stovalues {

double
my_f (const gsl_vector *v, void *params)
{
    auto args = static_cast<const Arguments *>(params);
    double sum = 0;

    for ( int i = 0 ; i < v->size ; i+=2 ) {
        sum += gsl_vector_get(v,i) * exp(gsl_vector_get(v, i+1));
    }

    logger()->trace("f called with N={}, sum is {}",v->size/2,sum);

    return sum;
}

/* The gradient of f, df = (df/dx, df/dy). */
void
my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
    auto args = static_cast<const Arguments *>(params);
    auto pi = std::numbers::pi_v<long double>;
    auto N = v->size/2;

    double x, y;
    double *p = (double *)params;

    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);

    gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
    gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
}

/* Compute both f and df together. */
void
my_fdf (const gsl_vector *x, void *params,
        double *f, gsl_vector *df)
{
    *f = my_f(x, params);
    my_df(x, params, df);
}


void Estimator::work(nlohmann::json &output_json)
{
    size_t iter = 0;
    int status;
    int N=1;
    int n = N*2;

    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;


    gsl_vector *x;
    gsl_multimin_function_fdf my_func;

    my_func.n = n;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = & ( this->args );

    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (n);
    gsl_vector_set (x, 0, 5.0);
    gsl_vector_set (x, 1, 7.0);

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
    while (status == GSL_CONTINUE && iter < 100);

    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);

}

} // namespace


