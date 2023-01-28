Include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.0.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

include (Catch)

include(CTest)

SET(CMAKE_CXX_FLAGS  "${CMAKE_C_FLAGS} ${GCC_WERR_COMPILE_FLAGS} -fsanitize=address -Werr")
SET(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(test-all  test-s2g.cpp ../src/environ.cpp ../src/logger.cpp ../src/parameters.cpp ../src/estimator.cpp ../src/beta_only_estimator.cpp)
find_package(spdlog REQUIRED)

target_include_directories(test-all PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries( test-all PUBLIC Catch2::Catch2WithMain spdlog::spdlog GSL::gsl GSL::gslcblas quadmath)


catch_discover_tests(test-all)
