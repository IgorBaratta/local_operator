#bin/bash

export CXX=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang++
export CC=/snx11273/home/ri-crichardson/spack/opt/spack/linux-rhel8-a64fx/gcc-12.2.0/llvm-16.0.0-hikjh2jya43tsrs4i35sduiswnccnl4p/bin/clang


export NUM_PROCS=48

run () {
    export CXX_FLAGS="-Ofast -march=armv8.2-a+sve -mcpu=a64fx -msve-vector-bits=512"
    for DEGREE in {1..8}
        do
            rm -rf build
            cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DBLOCK_SIZE=$3 .
            cmake --build build -j 10
        for i in {1..3}
        do
            mpirun -n ${NUM_PROCS} ./build/benchmark >> $4
            echo ", '${CXX_FLAGS}', ${NUM_PROCS}" >> $4
        done
    done
}

# Param :  1 - Precision 
#          2 - Batch Size 
#          3 - Block Size
#          4 - Output file name
run 8 1 0 "laplacian-tet-a64fx-llvm.txt"
run 8 1 1 "laplacian-tet-a64fx-llvm.txt"
run 8 1 16 "laplacian-tet-a64fx-llvm.txt"
run 8 8 0 "laplacian-tet-a64fx-llvm.txt"
run 8 8 1 "laplacian-tet-a64fx-llvm.txt"
run 8 8 16 "laplacian-tet-a64fx-llvm.txt"





# Param :  1 - Precision 
#          2 - Batch Size 
#          3 - Block Size
#          4 - Output file name
run 4 1 0 "laplacian-tet-a64fx-llvm.txt"
run 4 1 1 "laplacian-tet-a64fx-llvm.txt"
run 4 1 16 "laplacian-tet-a64fx-llvm.txt"
run 4 16 0 "laplacian-tet-a64fx-llvm.txt"
run 4 16 1 "laplacian-tet-a64fx-llvm.txt"
run 4 16 16 "laplacian-tet-a64fx-llvm.txt"
