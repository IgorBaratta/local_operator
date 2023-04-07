#bin/bash

export CXX=g++
export CXX_FLAGS="-march=native -Ofast"

export BATCH_SIZE=8
export BLOCK_SIZE=8
export PRECISION=8
export NUM_PROCS=1
export PRECOMPUTE=1

for DEGREE in 1
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} -DPRECOMPUTE=${PRECOMPUTE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> out.txt
    echo ", '${CXX_FLAGS}'" >> out.txt
done