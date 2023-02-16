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
    unsigned int number_of_terms, max_iterations, max_guesses;
    double beta_factor, initial_C, initial_beta;
    string guess_file="-";

    po::options_description desc("s2g parameters:");
    desc.add_options()
            ("number_of_terms,N", po::value<unsigned int>(&number_of_terms)->required(), "number of terms in the approximate sum")
            ("max_iterations,i", po::value<unsigned int>(&max_iterations)->default_value(204800), "Maximum number of iterations ")
            ("max_guesses,m", po::value<unsigned int>(&max_guesses)->default_value(0), "Maximum number of guesses")
            ("guess,g", po::value<string>(&guess_file), "If number_of_terms > 1, Initial Guess file from a previous run with N=N-1")
            ("beta_factor,f", po::value<double>(&beta_factor)->default_value(1.1), "beta growth factor")
            ("initial_C,C", po::value<double>(&initial_C)->default_value(1), "initial C")
            ("initial_beta,b", po::value<double>(&initial_beta)->default_value(0.1), "initial beta");

    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
        input_set["number_of_terms"] = number_of_terms;
        input_set["max_iterations"] = max_iterations;
        input_set["max_guesses"] = max_guesses;
        input_set["initial_C"] = initial_C;
        input_set["initial_beta"] = initial_beta;
        input_set["beta_factor"] = beta_factor;
        input_set["three_d"] = three_d;
        if ( not max_guesses and ( number_of_terms > 1 and vm.count("guess") < 1 ) ) {
            logger()->critical("FATAL: You must specify a guess file when number_of_terms > 1");
            exit(1);
        }
        input_set["guess"] = guess_file;

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

    //Calculated_C_Estimator constant_c_estimator(input_set);

    output_set["input"] = input_set;
    auto start_time = std::chrono::system_clock::now();

    program_info["start_time"] = time_to_string(start_time);

    Arguments args(input_set);

    if ( args.get_max_guesses() ) {
        double best_so_far = 1e20;
        nlohmann::json best_json;
        for (unsigned int i = 0U; i < args.get_max_guesses(); ++i) {
            Guess_Estimator estimator(args);
            estimator.minimize(output_set);

            if (estimator.get_estimate_error() < best_so_far) {
                best_so_far = estimator.get_estimate_error();
                best_json = estimator.get_output_set();
                cout << "Best :" << best_json << endl;
                output_set["best"]   = best_json;
            }
            output_set["best_method"] = "best";
        }
    } else {
        Incremental_Estimator estimator(args);
        estimator.minimize(output_set);
        output_set[estimator.get_method_name()] = estimator.get_output_set();
        output_set["best_method"] = estimator.get_method_name();
    }

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

    cout << output_set << endl;

    stovalues_global_cleanup();

}




