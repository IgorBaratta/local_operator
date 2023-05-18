#!/bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang

export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512 -DLIKWID_PERFMON"


export BATCH_SIZE=1
export PRECISION=8
export NUM_PROCS=48
export OPTIMIZE=1
export DEGREE=15


rm -rf output
mkdir -p output

cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .

cmake --build build -j 10

mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g MEM_DP -o output/MEM_DP_%h_%r.txt ./build/benchmark
mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g FLOPS_DP -o output/FLOPS_DP_%h_%r.txt ./build/benchmark
mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g L2 -o output/L2_%h_%r.txt ./build/benchmark
mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g DATA -o output/DATA_%h_%r.txt ./build/benchmark
mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g MEM -o output/MEM_%h_%r.txt ./build/benchmark
