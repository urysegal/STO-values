#include <boost/program_options.hpp>
#include <map>
#include <iostream>
#include <chrono>

#include "nlohmann/json.hpp"

#include "s2g.h"
#include "logger.h"
#include <sys/utsname.h>
#include <gsl/gsl_version.h>

namespace po = boost::program_options;
using namespace std;

namespace stovalues {

static nlohmann::json parse_args(int argc, const char **argv)
{
    nlohmann::json input_set;
    unsigned int number_of_terms, max_iterations, test_points;
    double accuracy, max_test_error;

    po::options_description desc("s2g parameters:");
    desc.add_options()
            ("number_of_terms,N", po::value<unsigned int>(&number_of_terms)->required(), "number of terms in the approximate sum")
            ("max_iterations", po::value<unsigned int>(&max_iterations)->default_value(1024), "Maximum number of iterations per trial N")
            ("accuracy", po::value<double>(&accuracy)->default_value(3.2e-8), "Accuracy desired")
            ("test_points", po::value<unsigned int>(&test_points)->default_value(1024), "Number of test points")
            ("max_test_error", po::value<double>(&max_test_error)->default_value(1e-10), "Maximum error allowed at any point");
    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        input_set["number_of_terms"] = number_of_terms;
        input_set["max_iterations"] = max_iterations;
        input_set["accuracy"] = accuracy;
        input_set["test_points"] = test_points;
        input_set["max_test_error"] = max_test_error;
    } catch (po::error &e) {
        std::stringstream desc_str;
        desc_str << desc;
        logger()->critical("FATAL: {}", std::string(e.what()) + "\n" + desc_str.str());
        exit(1);
    }

    return input_set;
}

static nlohmann::json parse_stdin()
{
    nlohmann::json input_set;
    try {
        cin >> input_set;
    } catch (exception &e) {
        logger()->critical("FATAL: {}", std::string(e.what()));
        exit(1);
    }
    return input_set;
}
}

using namespace stovalues;

static std::string time_to_string(const std::chrono::system_clock::time_point& tp)
{
    std::time_t t = std::chrono::system_clock::to_time_t(tp);
    return std::ctime(&t);
}

void
add_system_info(nlohmann::json &program_info)
{
    struct utsname name;
    uname(&name);
    program_info["num_cpus"] = sysconf( _SC_NPROCESSORS_ONLN );
    program_info["machine"] = string(name.machine);
    program_info["OS"] = string (name.sysname) + " " + string (name.release) + " " +string (name.version);
    program_info["compiler_version"] = __VERSION__;
    program_info["boost_version"] = std::to_string(BOOST_VERSION) ;
    program_info["gsl_version"] = GSL_VERSION;
}


void
estimate(nlohmann::json &input_set, nlohmann::json &output_set)
{
    nlohmann::json program_info;

    add_system_info(program_info);

    Estimator estimator(input_set);

    output_set["input"] = input_set;
    auto start_time = std::chrono::system_clock::now();

    program_info["start_time"] = time_to_string(start_time);

    estimator.minimize_C_beta_together(output_set);
    auto end_time = std::chrono::system_clock::now();

    program_info["end_time"] = time_to_string(end_time);
    program_info["run_time_seconds"] = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    output_set["program_info"] = program_info;
}


int
main(int argc, const char *argv[])
{
    stovalues_global_init();

    nlohmann::json input_set;
    nlohmann::json output_set;

    if ( argc > 1 ) {
        input_set = stovalues::parse_args(argc, argv);
    } else {
        input_set = stovalues::parse_stdin();
    }

    estimate(input_set, output_set);

    cout << output_set;

    stovalues_global_cleanup();

}




