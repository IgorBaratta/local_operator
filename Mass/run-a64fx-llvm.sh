#bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang

export NUM_PROCS=48

run () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
    for DEGREE in {1..8}
        do
            rm -rf build
            cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} .
            cmake --build build -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build/benchmark >> $3
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $3
        done
    done
}

# Order : Precision, Batch Size, Optimize Sum Factorization
run 8 8 "mass-tet-a64fx-llvm.txt"
run 8 1 "mass-tet-a64fx-llvm.txt"


# Single Precision
# Order : Precision, Batch Size, Optimize Sum Factorization
run 4 16 "mass-tet-a64fx-llvm.txt"
run 4 1 "mass-tet-a64fx-llvm.txt"
