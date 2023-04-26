#bin/bash

export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/gcc

export CXX_FLAGS="-march=native -Ofast -mprefer-vector-width=256"
export BATCH_SIZE=1
export BLOCK_SIZE=16
export PRECISION=8
export NUM_PROCS=1
export PRECOMPUTE=0
export OPTIMIZE=1
export DEGREE=5

rm -rf build
cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=${PRECISION} -DBATCH_SIZE=${BATCH_SIZE} -DDEGREE=${DEGREE} -DBLOCK_SIZE=${BLOCK_SIZE} -DPRECOMPUTE=${PRECOMPUTE} -DOPTIMIZE_SUM_FACTORIZATION=${OPTIMIZE} .
cmake --build build -j 10
mpirun -n ${NUM_PROCS} ./build/benchmark
echo -e "\n"