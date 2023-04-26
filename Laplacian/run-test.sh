#bin/bash

export CXX=g++
export CXX_FLAGS="-march=native -Ofast"

export NUM_PROCS=1
export BATCH_SIZE=1
export PRECISION=8
export BLOCK_SIZE=4
export DEGREE=5


rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark >> out.txt
echo ", '${CXX_FLAGS}'" >> out.txt