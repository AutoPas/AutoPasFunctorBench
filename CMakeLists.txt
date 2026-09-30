cmake_minimum_required(VERSION 3.14)
project(AP_Functor_Bench)

set(CMAKE_CXX_STANDARD 17)

# Default to Release build type for accurate benchmark timings
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
endif ()

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

target_include_directories(AP_Functor_Bench PRIVATE functor-variants)

# Pass autopas info to C++ code via a macro
target_compile_definitions(AP_Functor_Bench PRIVATE
        AUTOPAS_BRANCH="${AUTOPAS_BRANCH}"
        AUTOPAS_COMMIT="${AUTOPAS_COMMIT_HASH}"
)

# Optional Intel ITT API for Intel Advisor / VTune profiling
option(ENABLE_ITT "Enable Intel ITT API for selective profiling" OFF)
if (ENABLE_ITT)
    include(advisor)
endif ()