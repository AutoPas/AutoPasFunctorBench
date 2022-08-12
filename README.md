# AutoPas Functor Bench

Minimalistic program to analyze the functor kernel performance of AutoPas on different architectures.

## Requirements
* CMake
* git (for CMake FetchContent)
* C++17 compiler

## Compile Instructions

These are the compile instructions used on FUGAKU. Details might be subject to change according to your system.

### GCC
```bash
cmake                                                                                           \
    -DCMAKE_BUILD_TYPE=Release                                                                  \
    -DCMAKE_C_COMPILER=gcc                                                                      \
    -DCMAKE_CXX_COMPILER=g++                                                                    \
    -DCMAKE_C_FLAGS=-march=native                                                               \
    -DCMAKE_CXX_FLAGS=-march=native                                                             \
    -DAUTOPAS_FORMATTING_TARGETS=OFF                                                            \
    -DAUTOPAS_VECTOR_INSTRUCTIONS=NATIVE                                                        \
    ..

CC=gcc CXX=g++ make AutoPasFunctorBench -j12
```

### FCC
```bash
cmake                                                                                           \
    -DCMAKE_BUILD_TYPE=Release                                                                  \
    -DCMAKE_C_COMPILER=fcc                                                                      \
    -DCMAKE_CXX_COMPILER=FCC                                                                    \
    -DCMAKE_C_FLAGS="-Nclang -mcpu=a64fx -std=c17 -stdlib=libc++ -msve-vector-bits=512"         \
    -DCMAKE_CXX_FLAGS="-Nclang -mcpu=a64fx -std=c++17 -stdlib=libstdc++ -msve-vector-bits=512"  \
    -DAUTOPAS_FORMATTING_TARGETS=OFF                                                            \
    -DAUTOPAS_VECTOR_INSTRUCTIONS=DEFAULT                                                       \
    ..

CC='fcc -Nclang -mcpu=a64fx -std=c17 -stdlib=libc++'            \
CXX='FCC -Nclang -mcpu=a64fx -std=c++17 -stdlib=libc++'         \
make AutoPasFunctorBench -j12
```

## AutoPas Version
To select which version of AutoPas to benchmark change e.g. `GIT_TAG` in `cmake/modules/autopas.cmake`.
