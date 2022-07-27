# Local Finite ELement Operator Benchmarks

Depends on:
- FFCx: The FEniCSx Form Compiler
  -  https://github.com/FEniCS/ffcx.git

Optional dependencies:
- TSFC : The form compiler for the Firedrake project.
  - https://github.com/firedrakeproject/tsfc
- FFC: The FEniCS Form Compiler
  - https://bitbucket.org/fenics-project/ffc/

## Usage
```bash
usage: run.py [-h] [--form_compiler {ffcx,ffc,tsfc}] [--scalar_type {double,float}] [--problem {Laplacian,Mass,Elasticity,N1curl,Stokes}] [--conf CONF] [--degree DEGREE [DEGREE ...]] [--nrepeats NREPEATS]
              [--batch_size {None,4,8}] [--global_size GLOBAL_SIZE] [--action] [--mpi_size MPI_SIZE]

Run local assembly benchmark.

optional arguments:
  -h, --help            show this help message and exit
  --form_compiler {ffcx,ffc,tsfc}
                        Form Compiler to use (default: ffcx)
  --scalar_type {double,float}
                        Scalar type to use (default: double)
  --problem {Laplacian,Mass,Elasticity,N1curl,Stokes}
                        Problem to run (default: Laplacian)
  --conf CONF           Configuration file describing the compilers and flags. (default: compilers.yaml)
  --degree DEGREE [DEGREE ...]
                        Polynomial degree to evaluate the operators. (default: range(1, 4))
  --nrepeats NREPEATS   Number of times to run each experiment. (default: 3)
  --batch_size {None,4,8}
  --global_size GLOBAL_SIZE
                        Global number of dofs (assuming shared are dofs are duplicated). (default: 1000000.0)
  --action              Specify whether to run the problems with matrix free approach. (default: False)
  --mpi_size MPI_SIZE   The number of mpi processes to use. (default: 1)

```

## Compiler Configuration File
Example of compiler configuration file (compilers.yaml):
```yaml
gcc-11:
  version:
    - 11.1.0
  cpp:
    - /usr/bin/g++-11
  cc:
    - /usr/bin/gcc-11
  flags:
    - -Ofast -march=native -mprefer-vector-width=256

clang:
  version:
    - 12.0.1
  cpp:
    - /usr/bin/clang++
  cc:
    - /usr/bin/clang
  flags:
    - -Ofast -march=native -mprefer-vector-width=256
```

## Examples:
### Matrix free Weighted Laplacian degrees 5 and 6
```bash
python3 run.py --problem Lagrange  --degree 5 6 --form_compiler=ffcx --matrix_free
python3 run.py --problem Lagrange  --degree 5 6 --form_compiler=tsfc --matrix_free
python3 run.py --problem Lagrange  --degree 5 6 --form_compiler=ffc --matrix_free
```
### N1curl - curl curl - Rank 2
```bash
python3 run.py --problem N1curl --degree 5
```



## Data description

Results:
output/{Problem}.txt

```
table = [machine, kernel, compiler name, compiler flags, polyonomial degree, ffc opts, number of cells, local assemble time]
```
