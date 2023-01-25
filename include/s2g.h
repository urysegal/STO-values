#pragma once
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_erf.h>
#include "logger.h"
#include <numbers>
#include <quadmath.h>
#include "nlohmann/json.hpp"

namespace stovalues {


void stovalues_global_cleanup();
void stovalues_global_init();

class Arguments {
public:
    Arguments(const nlohmann::json &json);

    auto get_number_of_terms() const { return number_of_terms; }
    auto get_max_iterations() const { return max_iterations; }

private:
    unsigned int number_of_terms = 128;
    unsigned int max_iterations = 1024;
    unsigned int test_points = 1024;
    double accuracy = 3.2e-8;
    double max_test_error = 1e-10;
};

class Estimator {
public:

    Estimator(const Arguments &_args) : args(_args) {}
    void work(nlohmann::json &output_json);

    double average_error (const gsl_vector *v);
    void average_error_df (const gsl_vector *v, gsl_vector *df);

private:

    Arguments args;

    double diff_by_Ci(const gsl_vector *v, size_t i);
    double diff_by_bi(const gsl_vector *v, size_t i);

};

}