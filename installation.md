# Installation

## Installing ffcx and dependencies

```bash
python3 -m venv env/ffcx
source env/ffcx/bin/activate

python3 -m pip install git+https://github.com/FEniCS/ufl.git
python3 -m pip install git+https://github.com/FEniCS/basix.git
python3 -m pip install git+https://github.com/FEniCS/ffcx.git
python3 -m pip install pyyaml

```

## Installing ffcx and dependencies with spack

```
a
```

# Roofline with Intel Advisor

## Intel advisor commands

advisor --collect=survey --project-dir=./advi --search-dir src:r=. -- ./build/benchmark
advisor -collect tripcounts -flop -stacks --project-dir=./advi --search-dir src:r=. -- ./build/benchmark
advisor --report=roofline --with-stack --project-dir=./advi --report-output=./advi/out/roofline.html


## Intel advisor with MPI

mpirun -gtool "advisor --collect=survey --project-dir=./advi_results" -n 6 ./build/benchmark
mpirun -gtool "advisor -collect tripcounts -flop -stacks --search-dir src:r=. --project-dir=./advi:1-6"  -n 6 ./build/benchmark
advisor --report=roofline --with-stack --project-dir=./advi --report-output=./advi/out/roofline.html