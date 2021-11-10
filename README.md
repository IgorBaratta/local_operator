# Benchmark FFCx

Depends on:
- FFCx: The FEniCSx Form Compiler
  -  https://github.com/FEniCS/ffcx.git


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

## How to run:
### Lagrange - Stiffness Matrix
```bash
python3 run.py Lagrange
```
### N1curl - curl curl 
```bash
python3 run N1curl --degre 1 2 3
```

### Stokes Taylor-Hood mixed
```bash
python3 run_problem.py Stokes.ufl
```

Results:
{Name}.txt


## Add more compilers or flags:
Append compiler names and flags to lines 10 and 11 of [run_problem](https://github.com/IgorBaratta/benchmark_function/blob/main/run_problem.py)


## Plotting data:
```bash
python3 graph.py Stokes.txt 2
python3 graph.py Lagrange.txt 3
```
etc.

## Data description
```
table = [machine, kernel, compiler name, compiler flags, polyonomial degree, ffc opts, number of cells, local assemble time]
```
