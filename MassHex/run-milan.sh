#bin/bash

############################
# GCC
############################
export CXX=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/g++
export CC=/rds/user/ia397/hpc-work/spack/opt/spack/linux-rocky8-zen/gcc-8.4.0/gcc-12.2.0-vabb3h6ueyk5rgnkslaixjzbuoacjjxo/bin/gcc
export NUM_PROCS=128

run () {
    export CXX_FLAGS="-Ofast -march=native -mprefer-vector-width=256"
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
run 8 1 0 "mass-hex-milan-gcc.txt"
run 8 4 0 "mass-hex-milan-gcc.txt"
run 8 1 1 "mass-hex-milan-gcc-optimize.txt"
run 8 4 1 "mass-hex-milan-gcc-optimize.txt"



# Single Precision
# Order : Precision, Batch Size, Optimize Sum Factorization
run 4 1 0 "mass-hex-milan-gcc.txt"
run 4 8 0 "mass-hex-milan-gcc.txt"
run 4 1 1 "mass-hex-milan-gcc-optimize.txt"
run 4 8 1 "mass-hex-milan-gcc-optimize.txt"
