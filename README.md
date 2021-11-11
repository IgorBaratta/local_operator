# Benchmark Local Tensor Assembly

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
usage: run.py [-h] [--form_compiler {ffcx,ffc,tsfc}] [--problem {Lagrange,Elasticity,N1curl,Stokes}] [--conf CONF] [--degree DEGREE [DEGREE ...]] [--nrepeats {1,2,3,4,5,6,7,8,9,10}] [--matrix_free]

Run local assembly benchmark.

optional arguments:
  -h, --help            show this help message and exit
  --form_compiler {ffcx,ffc,tsfc}
                        Form Compiler to use (default: ffcx)
  --problem {Lagrange,Elasticity,N1curl,Stokes}
                        Problem to run (default: Lagrange)
  --conf CONF           Configuration file describing the compilers and flags. (default: compilers.yaml)
  --degree DEGREE [DEGREE ...]
                        Polynomial degree to evaluate the operators. (default: range(1, 4))
  --nrepeats {1,2,3,4,5,6,7,8,9,10}
                        Polynomial degree to evaluate the operators (default: 3)
  --matrix_free         Specify whether to run the problems with matrix freee approach. (default: False)

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
