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
            cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DOPTIMIZE_SUM_FACTORIZATION=$3 .
            cmake --build build -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build/benchmark >> $4
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $4
        done
    done
}

# Order : Precision, Batch Size, Optimize Sum Factorization
run 8 1 0 "out-icelake-gcc.txt"
run 8 8 0 "out-icelake-gcc.txt"
run 8 1 1 "out-icelake-gcc-optimize.txt"
run 8 8 1 "out-icelake-gcc-optimize.txt"



# Single Precision
# Order : Precision, Batch Size, Optimize Sum Factorization
run 4 1 0 "out-icelake-gcc.txt"
run 4 16 0 "out-icelake-gcc.txt"
run 4 1 1 "out-icelake-gcc-optimize.txt"
run 4 16 1 "out-icelake-gcc-optimize.txt"