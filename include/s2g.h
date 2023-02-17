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
    auto get_max_guesses() const { return max_guesses; }
    auto get_initial_beta() const { return  initial_beta; }
    auto get_initial_C() const { return  initial_C; }
    auto get_beta_factor() const { return beta_factor; }
    bool is_three_d() const { return three_d;}

private:
    unsigned int number_of_terms = 128;
    unsigned int max_iterations = 1024;
    std::string guess_file;
    unsigned int max_guesses = 0;
    double initial_beta = 0.1;
    double initial_C = 1;
    double beta_factor = 1.1;
    bool three_d = true;
};

struct result_term {
    double C;
    double beta;
};

class Guess_Estimator {
public:

    Guess_Estimator(const Arguments &_args) : args(_args) , N(args.get_number_of_terms()){}
    virtual ~Guess_Estimator();
    virtual void minimize(nlohmann::json &output_json);

    virtual double average_error (const gsl_vector *v);
    virtual void average_error_df (const gsl_vector *v, gsl_vector *df);
    virtual std::string get_method_name() const { return "C_conjugate_gradients"; }
    auto get_estimate_error() const { return estimate_error; }
    const nlohmann::json get_output_set() const { return output_set; }

protected:

    Arguments args;
    const unsigned int N=0;
    gsl_vector *x = nullptr;
    double estimate_error = 0;
    size_t iter = 0;
    std::vector<result_term> terms;
    gsl_multimin_fdfminimizer *s = nullptr;
    double step_size = 1e-1;
    double tolerance = 1e-2;
    double stop_gradient = 1e-12;
    nlohmann::json output_set;


    static real_t errfunc(real_t sqrt_beta);

    virtual double diff_by_Ci(const gsl_vector *v, size_t i);
    virtual double diff_by_bi(const gsl_vector *v, size_t i);

    virtual double GET_C(const gsl_vector *v, unsigned int i) const { return gsl_vector_get(v, (N + i)); }
    virtual double GET_beta(const gsl_vector *v, unsigned int i) const { return gsl_vector_get(v, (i)); }

    virtual void setup_initial_guess();
    void output_results(nlohmann::json &output_json, const gsl_vector *C_vector, const gsl_vector *beta_vector);

    virtual real_t calculate_adjuster (const gsl_vector *v) {return 1;}


};

class Three_D_Estimator : public Guess_Estimator {
public:
    Three_D_Estimator(const Arguments &_args) : Guess_Estimator(_args) {}

    virtual ~Three_D_Estimator();


    virtual double average_error (const gsl_vector *v) override;

protected:
    double diff_by_Ci(const gsl_vector *v, size_t i) override;
    double diff_by_bi(const gsl_vector *v, size_t i) override;
    real_t get_Bi(const gsl_vector *v, int i) const;
    real_t calculate_adjuster (const gsl_vector *v) override;

    void adjust_C(gsl_vector *v);


};

class Incremental_Estimator : public Guess_Estimator {
public:
    Incremental_Estimator(const Arguments &_args) : Guess_Estimator(_args) {}
    virtual ~Incremental_Estimator() ;

    virtual void minimize(nlohmann::json &output_json) override;

    double average_error (const gsl_vector *v) override;
    void average_error_df (const gsl_vector *v, gsl_vector *df) override;


protected:
    gsl_vector *adjusted_x = nullptr;
    void calculate_adjusted_x(const gsl_vector *v);


    void setup_initial_guess() override;

};

}