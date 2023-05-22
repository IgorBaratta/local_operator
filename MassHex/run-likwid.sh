#!/bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang

export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512 -DLIKWID_PERFMON"


export BATCH_SIZE=16
export PRECISION=4
export NUM_PROCS=48
export OPTIMIZE=1

mkdir -p likwid
for DEGREE in {1..15}
do
  rm -rf output
  mkdir -p output
  cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .
  cmake --build build -j 10
  mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g L2 -o output/L2_%h_%r.txt ./build/benchmark
  grep -e "L2 data volume" output/L2* >> likwid/mass_L2_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
  grep -e "L2 bandwidth" output/L2* >> likwid/mass_L2_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
  
  mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g MEM_DP -o output/MEM_DP_%h_%r.txt ./build/benchmark
  grep -e "Memory data volume" output/MEM_DP* >> likwid/mass_MEM_DP_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
  grep -e "Memory bandwidth" output/MEM_DP* >> likwid/mass_MEM_DP_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
  
  mpirun -n ${NUM_PROCS} likwid-perfctr -m -C E:N:1 -g FLOPS_DP -o output/FLOPS_DP_%h_%r.txt ./build/benchmark
  grep -e "FP_DP_FIXED_OPS_SPEC" output/FLOPS_DP* >> likwid/mass_FLOPS_DP_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
  grep -e "FP_DP_SCALE_OPS_SPEC" output/FLOPS_DP* >> likwid/mass_FLOPS_DP_${BATCH_SIZE}_${DEGREE}_${OPTIMIZE}_512.txt
done
