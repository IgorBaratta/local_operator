#bin/bash

############################
# llvm
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen3/gcc-12.2.0/llvm-16.0.0-eopnezxwmqi6mqkckzprcduet7sa6eqd/bin/clang++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen3/gcc-12.2.0/llvm-16.0.0-eopnezxwmqi6mqkckzprcduet7sa6eqd/bin/clang

export NUM_PROCS=128

run () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=256"
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


# Double Precision
# Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 8 1 1 16 0 "laplacian-hex-icelake-llvm.txt"
run 8 1 1 16 1 "laplacian-hex-icelake-llvm.txt"
run 8 1 1 1 1 "laplacian-hex-icelake-llvm.txt"
run 8 1 1 0 1 "laplacian-hex-icelake-llvm.txt"
run 8 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run 8 1 0 16 1 "laplacian-hex-icelake-llvm.txt"
run 8 4 1 1 0 "laplacian-hex-icelake-llvm.txt"
run 8 4 0 1 0 "laplacian-hex-icelake-llvm.txt"


# Single Precision
# Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 4 1 1 16 0 "laplacian-hex-icelake-llvm.txt"
run 4 1 1 16 1 "laplacian-hex-icelake-llvm.txt"
run 4 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run 4 1 0 16 1 "laplacian-hex-icelake-llvm.txt"
run 4 8 1 1 0 "laplacian-hex-icelake-llvm.txt"
run 4 8 0 1 0 "laplacian-hex-icelake-llvm.txt"

