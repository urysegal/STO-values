#pragma once

#include "nlohmann/json.hpp"

namespace stovalues {


void stovalues_global_cleanup();
void stovalues_global_init();

class Arguments {
public:
    Arguments(const nlohmann::json &json);

private:
    unsigned int max_number_of_terms = 128;
    unsigned int max_iterations = 1024;
    unsigned int test_points = 1024;
    double accuracy = 3.2e-8;
    double max_test_error = 1e-10;
};

class Estimator {
public:

    Estimator(const Arguments &_args) : args(_args) {}
    void work(nlohmann::json &output_json);

private:

    Arguments args;

};

}