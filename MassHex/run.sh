#bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-10.3.0/gcc-12.2.0-zf65tfzmboe3t6lkegf4rrw2zwnprvtn/bin/g++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-10.3.0/gcc-12.2.0-zf65tfzmboe3t6lkegf4rrw2zwnprvtn/bin/gcc


# export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512 -mcmodel=large -fopt-info-vec"
export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=256"
export BATCH_SIZE=1
export PRECISION=8
export NUM_PROCS=48
export OPTIMIZE=0
export DEGREE=5

rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark

# mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L2 -g L2CACHE -g FLOPS_DP -o output/output256_%h_%r.txt ./build/benchmark


# mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L2 -o output/L2_%h_%r.txt -m ./build/benchmark
# echo -e "\n"
# mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L3 -o output/L3_%h_%r.txt -m ./build/benchmark
# echo -e "\n"
# mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g CLOCK -o output/CLOCK_%h_%r.txt -m ./build/benchmark
# echo -e "\n"




