#bin/bash

export CXX=g++
export CXX_FLAGS="-march=native -Ofast"
export BATCH_SIZE=4
export PRECISION=8
export NUM_PROCS=1
export DEGREE=2

rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark 1