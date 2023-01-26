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
    std::string guess_file;
};

class Estimator {
public:

    Estimator(const Arguments &_args) : args(_args) , N(args.get_number_of_terms()){}
    virtual void minimize(nlohmann::json &output_json);

    double average_error (const gsl_vector *v);
    void average_error_df (const gsl_vector *v, gsl_vector *df);

protected:

    Arguments args;
    const unsigned int N=0;
    gsl_vector *x = nullptr;

    double diff_by_Ci(const gsl_vector *v, size_t i);
    double diff_by_bi(const gsl_vector *v, size_t i);

    virtual double GET_C(const gsl_vector *v, unsigned int i) const { return gsl_vector_get(v, (N + i)); }

    void update_C(gsl_vector *);

};

class Calculated_C_Estimator : public Estimator
{
public:
    Calculated_C_Estimator(const Arguments &_args) : Estimator(_args) {}
    void average_error_df_beta_only(const gsl_vector *v, gsl_vector *df);
    void minimize(nlohmann::json &output_json) override;

protected:
    double GET_C(const gsl_vector *v, unsigned int i) const override { return gsl_vector_get(x, (N + i)); }

};

}