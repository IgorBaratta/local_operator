#bin/bash

############################
# llvm
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-a64fx/llvm-11.2.0/llvm-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-a64fx/llvm-11.2.0/llvm-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/llvm

export NUM_PROCS=48

run512 () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512 -mcmodel=large -fopt-info-vec"
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

run256 () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=256"
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
run512 8 1 0 "mass-hex-a64fx-llvm.txt"
run512 8 8 0 "mass-hex-a64fx-llvm.txt"
run512 8 1 1 "mass-hex-a64fx-llvm.txt"
run256 8 1 1 "mass-hex-a64fx-llvm.txt"

# Single Precision
run512 4 1 0 "mass-hex-a64fx-llvm.txt"
run512 4 16 0 "mass-hex-a64fx-llvm.txt"
run512 4 1 1 "mass-hex-a64fx-llvm.txt"
run256 4 1 1 "mass-hex-a64fx-llvm.txt"
