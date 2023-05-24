#bin/bash

export CXX=g++
export CC=gcc

export BATCH_SIZE=16
export PRECISION=8
export NUM_PROCS=1
export DEGREE=5
export CXX_FLAGS="-march=native -Ofast -mprefer-vector-width=256 -DNDEBUG"
export BLIS_DIR=/home/igorbaratta/Projects/spack/opt/spack/linux-ubuntu22.04-skylake/gcc-11.3.0/blis-0.9.0-upogwxhzgzmw2oltilp2yakehonkspkk


rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark