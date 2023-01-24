#include <iostream>
#include "nlohmann/json.hpp"
#include "s2g.h"
#include "logger.h"

using namespace nlohmann;

namespace stovalues
{

Arguments::Arguments(const nlohmann::json &j)
{
    try {
        this->accuracy = j["accuracy"];
        this->max_iterations = j["max_iterations"];
        this->max_number_of_terms = j["max_number_of_terms"];
        this->test_points = j["test_points"];
        this->max_test_error = j["max_test_error"];
    } catch (json::parse_error& ex) {
        std::cerr << "parse error : " << ex.what() << std::endl;
    }
}


}
