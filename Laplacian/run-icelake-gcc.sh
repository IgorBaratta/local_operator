#bin/bash

############################
# GCC
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-icelake/gcc-11.2.0/gcc-12.2.0-ipuhem4nelxd7n2us6tcshvwzkce5ip6/bin/gcc
export NUM_PROCS=76

run512 () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=512"
    for DEGREE in {1..8}
        do
            rm -rf build
            cmake -B build/ -DCMAKE_CXX_FLAGS="${CXX_FLAGS}" -DPRECISION=$1 -DBATCH_SIZE=$2 -DDEGREE=${DEGREE} -DBLOCK_SIZE=$3.
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
run512 8 1 0 "laplacian-tet-icelake-gcc.txt"
run512 8 1 1 "laplacian-tet-icelake-gcc.txt"
run512 8 1 16 "laplacian-tet-icelake-gcc.txt"
run512 8 8 0 "laplacian-tet-icelake-gcc.txt"
run512 8 8 1 "laplacian-tet-icelake-gcc.txt"
run512 8 8 16 "laplacian-tet-icelake-gcc.txt"

# Param :  1 - Precision 
#          2 - Batch Size 
#          3 - Block Size
#          4 - Output file name
run512 4 1 0 "laplacian-tet-icelake-gcc.txt"
run512 4 1 1 "laplacian-tet-icelake-gcc.txt"
run512 4 1 16 "laplacian-tet-icelake-gcc.txt"
run512 4 16 0 "laplacian-tet-icelake-gcc.txt"
run512 4 16 1 "laplacian-tet-icelake-gcc.txt"
run512 4 16 16 "laplacian-tet-icelake-gcc.txt"
