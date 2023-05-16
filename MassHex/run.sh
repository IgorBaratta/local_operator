#bin/bash

export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/gcc

export CXX_FLAGS="-march=native -Ofast -mprefer-vector-width=256 -DLIKWID_PERFMON"
export BATCH_SIZE=1
export PRECISION=8
export NUM_PROCS=76
export OPTIMIZE=1
export DEGREE=15

rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .
cmake --build build -j 10
# mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L2 -g L2CACHE -g FLOPS_DP -o output/output256_%h_%r.txt ./build/benchmark
mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L2 -o output/L2_%h_%r.txt -m ./build/benchmark
echo -e "\n"
mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g L3 -o output/L3_%h_%r.txt -m ./build/benchmark
echo -e "\n"
mpirun -n ${NUM_PROCS} likwid-perfctr -C E:N:1 -g CLOCK -o output/CLOCK_%h_%r.txt -m ./build/benchmark
echo -e "\n"




