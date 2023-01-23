#include <memory>
#include <iostream>
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "logger.h"

namespace stovalues {

class Logger {
private:
    static constexpr const char *log_file_name = "s2g.log";
    std::shared_ptr<spdlog::logger> logger;
    std::shared_ptr<spdlog::sinks::stderr_color_sink_mt> console_sink;
    std::shared_ptr<spdlog::sinks::basic_file_sink_mt> logfile_sink;
public:

    spdlog::logger *get_logger() { return logger.get();}

    void setup_logger() {

        try {
            console_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
            console_sink->set_level(spdlog::level::err);
            console_sink->set_pattern("[%^%l%$] %v");

            const char *LOG_DIR = getenv("LOG_DIR");
            std::string logfile = std::string(LOG_DIR? LOG_DIR : "/tmp" ) + "/" + log_file_name;
            logfile_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logfile, true);
            logfile_sink->set_level(spdlog::level::trace);

            auto logger_ = new spdlog::logger("s2g", {console_sink, logfile_sink});
            logger.reset(logger_);
            logger->set_level(spdlog::level::trace);
            logger->info("Logger started");
        } catch ( std::exception &e){
            std::cerr << "Cannot set up logger: " << e.what(); // LCOV_EXCL_LINE
            exit(1); // LCOV_EXCL_LINE
        } // LCOV_EXCL_LINE
    }
};

static Logger *global_logger = nullptr;
static std::mutex logger_create_lock;

spdlog::logger *logger()
{
    std::lock_guard<std::mutex> guard(logger_create_lock);
    if ( not global_logger )
    {
        global_logger = new Logger();
        global_logger->setup_logger();
    }
    return global_logger->get_logger();
}

}
