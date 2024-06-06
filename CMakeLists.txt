cmake_minimum_required(VERSION 3.14)
project(AutoPasFunctorBench)

set(CMAKE_CXX_STANDARD 17)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)
include(autopas)

# TODO fcc enables this by default. Add flags to (de-)activate depending on the compiler
option(ENABLE_FAPP "Activate calls to Fujitsu Application Performance Profiler" OFF)

option(FORCE_AVX2 "Forcibly downgrades functor to AVX2, even if on machine with AVX512 enabled" OFF)

option(ENABLE_HWY "Activates HWY functor" ON)

option(ENABLE_XSIMD "Activates XSIMD functor" ON)

option(ENABLE_SIMDe "Activates SIMDe functor" ON)

option(ENABLE_MIPP "Activates MIPP functor" ON)

option(ENABLE_AUTOVEC "Activates Autovec functor" ON)

add_executable(AutoPasFunctorBench main.cpp)
target_link_libraries(AutoPasFunctorBench
        PUBLIC
        autopas
        molecularDynamicsLibrary
        )
target_compile_definitions(AutoPasFunctorBench
        PUBLIC
        $<$<BOOL:${ENABLE_FAPP}>:ENABLE_FAPP>
        $<$<BOOL:${FORCE_AVX2}>:FORCE_AVX2>
        $<$<BOOL:${ENABLE_HWY}>:ENABLE_HWY>
        $<$<BOOL:${ENABLE_XSIMD}>:ENABLE_XSIMD>
        $<$<BOOL:${ENABLE_SIMDe}>:ENABLE_SIMDe>
        $<$<BOOL:${ENABLE_MIPP}>:ENABLE_MIPP>
        $<$<BOOL:${ENABLE_AUTOVEC}>:ENABLE_AUTOVEC>
        )