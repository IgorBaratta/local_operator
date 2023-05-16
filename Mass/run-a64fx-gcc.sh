#bin/bash

export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/gcc

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
run 8 4 "mass-tet-a64fx-gcc.txt"
run 8 1 "mass-tet-a64fx-gcc.txt"


# Single Precision
# Order : Precision, Batch Size, Optimize Sum Factorization
run 4 1 "mass-tet-a64fx-gcc.txt"
run 4 8 "mass-tet-a64fx-gcc.txt"
