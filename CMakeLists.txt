cmake_minimum_required(VERSION 3.14)
project(AP_Functor_Bench)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
include(google-benchmark)
include(autopas)
include(cli11)

add_executable(AP_Functor_Bench main.cpp)
target_link_libraries(AP_Functor_Bench
        PRIVATE
        benchmark::benchmark
        CLI11::CLI11
        PUBLIC
        autopas
        molecularDynamicsLibrary
        )

# Pass autopas info to C++ code via a macro
target_compile_definitions(AP_Functor_Bench PRIVATE
        AUTOPAS_BRANCH="${AUTOPAS_BRANCH}"
        AUTOPAS_COMMIT="${AUTOPAS_COMMIT_HASH}"
)