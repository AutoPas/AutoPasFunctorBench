cmake_minimum_required(VERSION 3.14)
project(AutoPasFunctorValidation)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
include(autopas)

add_executable(AutoPasFunctorValidation main.cpp)
target_link_libraries(AutoPasFunctorValidation
        PUBLIC
        autopas
        molecularDynamicsLibrary
        )