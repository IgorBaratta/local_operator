#bin/bash

############################
# GCC
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/gcc


export NUM_PROCS=128

run () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=256"
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
run 8 1 0 "laplacian-tet-milan-gcc.txt"
run 8 1 1 "laplacian-tet-milan-gcc.txt"
run 8 1 16 "laplacian-tet-milan-gcc.txt"
run 8 4 0 "laplacian-tet-milan-gcc.txt"
run 8 4 1 "laplacian-tet-milan-gcc.txt"
run 8 4 16 "laplacian-tet-milan-gcc.txt"

# Param :  1 - Precision 
#          2 - Batch Size 
#          3 - Block Size
#          4 - Output file name
run 4 1 0 "laplacian-tet-milan-gcc.txt"
run 4 1 1 "laplacian-tet-milan-gcc.txt"
run 4 1 16 "laplacian-tet-milan-gcc.txt"
run 4 8 0 "laplacian-tet-milan-gcc.txt"
run 4 8 1 "laplacian-tet-milan-gcc.txt"
run 4 8 16 "laplacian-tet-milan-gcc.txt"