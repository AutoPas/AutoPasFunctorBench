# autopas library
message(STATUS "Adding Google Benchmark.")

# Enable ExternalProject CMake module
include(FetchContent)

set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_USE_BUNDLED_GTEST OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY "https://github.com/google/benchmark.git"
        GIT_TAG "v1.9.4"
)
# Populate dependency
FetchContent_MakeAvailable(googlebenchmark)

#target_link_libraries(AutoPasFunctorValidation PRIVATE benchmark::benchmark)