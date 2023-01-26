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
        this->max_iterations = j["max_iterations"];
        this->number_of_terms = j["number_of_terms"];
        this->guess_file = j["guess"];
    } catch (json::parse_error& ex) {
        std::cerr << "parse error : " << ex.what() << std::endl;
    }
}


}
