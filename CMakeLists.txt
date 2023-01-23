cmake_minimum_required(VERSION 3.16)

project(stovalues)

include(CTest)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

option(COVERAGE "Coverage" ON)

set(COMPILE_FLAGS -Wall -Werror -fsanitize=address -fsanitize=leak -fno-omit-frame-pointer -pthread)
set(LINK_FLAGS -fsanitize=address)


find_package(Boost 1.71 REQUIRED program_options)


find_package(PkgConfig)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(${stovalues_SOURCE_DIR} ${stovalues_SOURCE_DIR}/include)



find_package(spdlog REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(s2g src/s2g.cpp src/environ.cpp src/logger.cpp)
add_dependencies(s2g nlohmann_json)
target_link_libraries(s2g PUBLIC spdlog::spdlog)
target_link_libraries (s2g PUBLIC Eigen3::Eigen)
target_link_libraries (s2g PUBLIC ${Boost_LIBRARIES})
target_include_directories(s2g PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(s2g PUBLIC ${COMPILE_FLAGS})
target_link_options(s2g PUBLIC ${LINK_FLAGS})

add_subdirectory(test)
