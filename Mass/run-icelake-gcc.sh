#bin/bash

############################
# GCC
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/gcc

run512 () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=512"
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

# Order : Precision, Batch Size, Name
run512 8 1 "mass-tet-icelake-gcc.txt"
run512 8 8 "mass-tet-icelake-gcc.txt"

# Single Precision
# Order : Precision, Batch Size, Name
run512 4 1 "mass-tet-icelake-gcc.txt"
run512 4 16 "mass-tet-icelake-gcc.txt"