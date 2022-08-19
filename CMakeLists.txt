cmake_minimum_required(VERSION 3.14)
project(AutoPasFunctorBench)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
include(autopas)

# TODO fcc enables this by default. Add flags to (de-)activate depending on the compiler
option(ENABLE_FAPP "Activate calls to Fujitsu Application Performance Profiler" OFF)

add_executable(AutoPasFunctorBench main.cpp)
target_link_libraries(AutoPasFunctorBench
        PUBLIC
        autopas
        )
target_compile_definitions(AutoPasFunctorBench
        PUBLIC
        $<$<BOOL:${ENABLE_FAPP}>:ENABLE_FAPP>
        )