#!/bin/bash

# Simulation parameters
PARTICLES=5000
ITERATIONS=2
DENSITY=0.5
OUTPUTFOLDERNAME=measurePerf_

# Create output directory
mkdir -p ${OUTPUTFOLDERNAME}

# Loop over precisions and types
for PRECISION in SPSP SPDP DPDP; do
	mkdir -p build${PRECISION}
	cd build${PRECISION}
	for TYPE in single pair; do
        echo "Running ${TYPE} with precision ${PRECISION}"

		CC=clang CXX=clang++ cmake ..
        
        # Configure CMake with current precision
        if ! cmake -DAUTOPAS_PRECISION_MODE=${PRECISION} -DAUTOPAS_VECTOR_INSTRUCTIONS=AVX2 .; then
            echo "CMake configuration failed for ${PRECISION}"
            continue
        fi
        
        # Build the project
        if ! cmake --build . --parallel 8; then
            echo "Build failed for ${PRECISION}"
            continue
        fi

        #give CPU some time to breath
        echo "--------------------sleeping for 180--------------------"
        sleep 180
        echo "---------------------slept for 180---------------------"
        
        # Run benchmark and generate output file
        if ! ./AutoPasFunctorBench "${TYPE}" "${ITERATIONS}" "${PARTICLES}" "${PARTICLES}" "${DENSITY}" "../${OUTPUTFOLDERNAME}/${TYPE}${PRECISION}.json"; then
            echo "Benchmark failed for ${TYPE} with ${PRECISION}"
            continue
        fi
        
        echo "Completed ${TYPE} with precision ${PRECISION}"
        echo "----------------------------------------"
    done
	cd ..
done

echo "All benchmarks completed. Results are in the ${OUTPUTFOLDERNAME} directory."