#bin/bash

############################
# GCC
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/gcc
export NUM_PROCS=76

export PRECISION=8
export BATCH_SIZE=1

run () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=512"
    for DEGREE in {1..15}
        do
            rm -rf build
            cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DPRECOMPUTE=$3 -DBLOCK_SIZE=$4 -DOPTIMIZE_SUM_FACTORIZATION=$5 .
            cmake --build build -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build/benchmark >> $6
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $6
        done
    done
}

# Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 8 1 1 16 0 "out-icelake-gcc-precompute.txt"
run 8 1 1 16 1 "out-icelake-gcc-precompute-optimize.txt"

run 8 1 0 16 0 "out-icelake-gcc-no-precompute.txt"
run 8 1 0 16 1 "out-icelake-gcc-no-precompute-optimize.txt"

run 8 8 1 1 0 "out-icelake-gcc-precompute.txt"
run 8 8 0 1 0 "out-icelake-gcc-no-precompute.txt"


