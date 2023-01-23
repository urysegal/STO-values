#include "s2g.h"
#include "logger.h"
#include <boost/program_options.hpp>
#include <map>
#include "nlohmann/json.hpp"



namespace po = boost::program_options;
using namespace std;

namespace stovalues {

static nlohmann::json parse_args(int argc, const char **argv)
{
    nlohmann::json input_set;

    std::string option_basis_set, option_xyz_URL;

    po::options_description desc("System information parameters:");
    desc.add_options()
            ("xyz", po::value<std::string>(&option_xyz_URL)->required(), "URL of an XYZ file of the system")
            ("basis_set", po::value<std::string>(&option_basis_set)->required(), "Basis Set");
    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc - 1, &argv[1], desc), vm);
        po::notify(vm);


    } catch (po::error &e) {
        std::stringstream desc_str;
        desc_str << desc;
        logger()->critical("FATAL: {}", std::string(e.what()) + "\n" + desc_str.str());
    }

    return input_set;
}

}


int
main(int argc, const char *argv[])
{
    stovalues::stovalues_global_cleanup();
    nlohmann::json input_set;

    if ( argc == 1 ) {
        input_set = stovalues::parse_args(argc, argv);
    } else {

    }

    stovalues::stovalues_global_init();

}




