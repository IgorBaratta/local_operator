#bin/bash

export CXX=clang++-15
export CXX_FLAGS="-march=native -Ofast"

export BATCH_SIZE=4
export BLOCK_SIZE=1
export PRECISION=8
export NUM_PROCS=1
export PRECOMPUTE=1
export OPTIMIZE=1
export DEGREE=5

rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} -DPRECOMPUTE=${PRECOMPUTE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark >> out.txt
echo ", '${CXX_FLAGS}'" >> out.txt
