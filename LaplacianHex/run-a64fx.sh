#bin/bash

export CXX=g++
export CC=gcc
export CXX_FLAGS="-Ofast -march=native -mcpu=a64fx -mcmodel=large -fopt-info-vec -msve-vector-bits=512"
export OUTPUT="out-a64fx-gcc.txt"

export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=8
export BLOCK_SIZE=8

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=8
export PRECISION=8
export BLOCK_SIZE=1

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15	      
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=4
export BLOCK_SIZE=16

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=16
export PRECISION=4
export BLOCK_SIZE=1

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15	      
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

# repeat for clang

export CXX=clang++
export CC=clang
export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
export OUTPUT="out-a64fx-clang.txt"

export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=8
export BLOCK_SIZE=8

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=8
export PRECISION=8
export BLOCK_SIZE=1

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15	      
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=1
export PRECISION=4
export BLOCK_SIZE=16

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done

export NUM_PROCS=48
export BATCH_SIZE=16
export PRECISION=4
export BLOCK_SIZE=1

for DEGREE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15	      
do
    rm -rf build
    cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} .
    cmake --build build -j 10
    mpirun -n ${NUM_PROCS} ./build/benchmark >> ${OUTPUT}
    echo ", '${CXX_FLAGS}'" >> ${OUTPUT}
done
