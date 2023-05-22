#bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang
export NUM_PROCS=48


run () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
    for DEGREE in {1..15}
        do
            rm -rf build_a64fx
            cmake -B build_a64fx/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DPRECOMPUTE=$3 -DBLOCK_SIZE=$4 -DOPTIMIZE_SUM_FACTORIZATION=$5 .
            cmake --build build_a64fx -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build_a64fx/benchmark >> $6
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $6
        done
    done
}


# Double Precision
# Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 8 1 1 16 0 "laplacian-hex-a64fx-llvm.txt"
run 8 1 1 16 1 "laplacian-hex-a64fx-llvm.txt"
run 8 1 0 16 0 "laplacian-hex-a64fx-llvm.txt"
run 8 1 0 16 1 "laplacian-hex-a64fx-llvm.txt"
run 8 8 1 1 0 "laplacian-hex-a64fx-llvm.txt"
run 8 8 0 1 0 "laplacian-hex-a64fx-llvm.txt"


# Single Precision
## Order : Precision, Batch Size, PRECOMPUTE, Block Size, Optimize Sum Factorization
run 4 1 1 16 0 "laplacian-hex-a64fx-llvm.txt"
run 4 1 1 16 1 "laplacian-hex-a64fx-llvm.txt"
run 4 1 0 16 0 "laplacian-hex-a64fx-llvm.txt"
run 4 1 0 16 1 "laplacian-hex-a64fx-llvm.txt"
run 4 16 1 1 0 "laplacian-hex-a64fx-llvm.txt"
run 4 16 0 1 0 "laplacian-hex-a64fx-llvm.txt"