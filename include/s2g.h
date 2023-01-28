#pragma once
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>
#include "logger.h"
#include <numbers>
#include <quadmath.h>
#include "nlohmann/json.hpp"

namespace stovalues {

typedef  __float128 real_t;


void stovalues_global_cleanup();
void stovalues_global_init();

class Arguments {
public:
    Arguments(const nlohmann::json &json);

    auto get_number_of_terms() const { return number_of_terms; }
    auto get_max_iterations() const { return max_iterations; }
    auto get_guess_file() const { return guess_file; }
private:
    unsigned int number_of_terms = 128;
    unsigned int max_iterations = 1024;
    std::string guess_file;
};

struct result_term {
    double C;
    double beta;
};

class Estimator {
public:

    Estimator(const Arguments &_args) : args(_args) , N(args.get_number_of_terms()){}
    virtual ~Estimator();
    virtual void minimize(nlohmann::json &output_json);

    double average_error (const gsl_vector *v);
    void average_error_df (const gsl_vector *v, gsl_vector *df);
    virtual std::string get_method_name() const { return "C_conjugate_gradients"; }
    auto get_estimate_error() const { return estimate_error; }

protected:

    Arguments args;
    const unsigned int N=0;
    gsl_vector *x = nullptr;
    double estimate_error = 0;
    size_t iter = 0;
    std::vector<result_term> terms;
    gsl_multimin_fdfminimizer *s = nullptr;
    double step_size = 1e-6;
    double tolerance = 1e-11;
    double stop_gradient = 1e-12;

    static real_t errfunc(real_t sqrt_beta);


    double diff_by_Ci(const gsl_vector *v, size_t i);
    double diff_by_bi(const gsl_vector *v, size_t i);

    virtual double GET_C(const gsl_vector *v, unsigned int i) const { return gsl_vector_get(v, (N + i)); }

    void update_C(gsl_vector *);
    void setup_initial_guess();
    void output_results(nlohmann::json &output_json, const gsl_vector *C_vector, const gsl_vector *beta_vector);


};

class Calculated_C_Estimator : public Estimator
{
public:
    Calculated_C_Estimator(const Arguments &_args) : Estimator(_args) {}
    void average_error_df_beta_only(const gsl_vector *v, gsl_vector *df);
    void minimize(nlohmann::json &output_json) override;
    std::string get_method_name() const override { return "C_implied"; }

protected:
    double GET_C(const gsl_vector *v, unsigned int i) const override { return gsl_vector_get(x, (N + i)); }

};

}