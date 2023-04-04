#bin/bash

export CXX=g++
export CXX_FLAGS="-Ofast -march=native -mcpu=a64fx -mcmodel=large -fopt-info-vec -msve-vector-bits=512"

export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=8
export BLOCK_SIZE=8

for DEGREE in 1 2 3 4 5 6 7 8
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> out-a64fx.txt
    echo ", '${CXX_FLAGS}'" >> out-a64fx.txt
done

export NUM_PROCS=48
export BATCH_SIZE=8
export PRECISION=8
export BLOCK_SIZE=1

for DEGREE in 1 2 3 4 5 6 7 8
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> out-a64fx.txt
    echo ", '${CXX_FLAGS}'" >> out-a64fx.txt
done


export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=4
export BLOCK_SIZE=16

for DEGREE in 1 2 3 4 5 6 7 8
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> out-a64fx.txt
    echo ", '${CXX_FLAGS}'" >> out-a64fx.txt
done

export NUM_PROCS=48
export BATCH_SIZE=16
export PRECISION=4
export BLOCK_SIZE=1

for DEGREE in 1 2 3 4 5 6 7 8
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> out-a64fx.txt
    echo ", '${CXX_FLAGS}'" >> out-a64fx.txt
done
