#bin/bash

############################
# LLVM
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/llvm-16.0.0-7o7njhjb7ggirmrmnddky3ymb2vopwn3/bin/clang++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/llvm-16.0.0-7o7njhjb7ggirmrmnddky3ymb2vopwn3/bin/clang

export NUM_PROCS=76

run () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=512"
    for DEGREE in {1..15}
        do
            rm -rf build_icelake
            cmake -B build_icelake/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DPRECOMPUTE=$3 -DBLOCK_SIZE=$4 -DOPTIMIZE_SUM_FACTORIZATION=$5 .
            cmake --build build_icelake -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build_icelake/benchmark >> $6
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $6
        done
    done
}


run256 () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=256"
    for DEGREE in {1..15}
        do
            rm -rf build_icelake
            cmake -B build_icelake/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DPRECOMPUTE=$3 -DBLOCK_SIZE=$4 -DOPTIMIZE_SUM_FACTORIZATION=$5 .
            cmake --build build_icelake -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build_icelake/benchmark >> $6
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $6
        done
    done
}

# Double Precision
# Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 8 1 1 16 0 "laplacian-hex-icelake-llvm.txt"
run 8 1 1 16 1 "laplacian-hex-icelake-llvm.txt"
run 8 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run 8 1 0 16 1 "laplacian-hex-icelake-llvm.txt"
run 8 8 1 1 0 "laplacian-hex-icelake-llvm.txt"
run 8 8 0 1 0 "laplacian-hex-icelake-llvm.txt"


# Single Precision
## Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 4 1 1 16 0 "laplacian-hex-icelake-llvm.txt"
run 4 1 1 16 1 "laplacian-hex-icelake-llvm.txt"
run 4 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run 4 1 0 16 1 "laplacian-hex-icelake-llvm.txt"
run 4 16 1 1 0 "laplacian-hex-icelake-llvm.txt"
run 4 16 0 1 0 "laplacian-hex-icelake-llvm.txt"


## Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run256 8 1 1 1 0 "laplacian-hex-icelake-llvm.txt"
run256 8 1 1 1 1 "laplacian-hex-icelake-llvm.txt"
run256 8 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run256 8 1 0 16 1 "laplacian-hex-icelake-llvm.txt"

## Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run256 4 1 1 16 0 "laplacian-hex-icelake-llvm.txt"
run256 4 1 1 16 1 "laplacian-hex-icelake-llvm.txt"
run256 4 1 0 16 0 "laplacian-hex-icelake-llvm.txt"
run256 4 1 0 16 1 "laplacian-hex-icelake-llvm.txt"