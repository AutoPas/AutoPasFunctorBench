cmake_minimum_required(VERSION 3.14)
project(AutoPasFunctorValidation)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
include(google-benchmark)
include(autopas)

add_executable(AutoPasFunctorValidation main.cpp)
target_link_libraries(AutoPasFunctorValidation
        PRIVATE
        benchmark::benchmark
        PUBLIC
        autopas
        molecularDynamicsLibrary
        )

# Pass autopas info to C++ code via a macro
target_compile_definitions(AutoPasFunctorValidation PRIVATE
        AUTOPAS_BRANCH="${AUTOPAS_BRANCH}"
        AUTOPAS_COMMIT="${AUTOPAS_COMMIT_HASH}"
)