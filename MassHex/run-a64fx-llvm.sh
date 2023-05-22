#bin/bash

############################
# llvm
############################
export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang

export NUM_PROCS=48

run512 () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
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
# run512 8 1 0 "mass-hex-a64fx-llvm.txt"
# run512 8 8 0 "mass-hex-a64fx-llvm-8.txt"
# run512 8 1 1 "mass-hex-a64fx-llvm.txt"
# run256 8 1 1 "mass-hex-a64fx-llvm.txt"

# Single Precision
# run512 4 1 0 "mass-hex-a64fx-llvm.txt"
run512 4 16 0 "mass-hex-a64fx-llvm-16.txt"
# run512 4 1 1 "mass-hex-a64fx-llvm.txt"
# run256 4 1 1 "mass-hex-a64fx-llvm.txt"
