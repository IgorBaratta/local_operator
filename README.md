# Benchmark FFCx

Depends on:
- ffcx -  branch https://github.com/FEniCS/ffcx/tree/igor/loop-invariant


## How to run:
### Lagrange - Stiffness Matrix
```bash
python3 run_problem.py Lagrange.ufl
```
### N1curl - curl curl 
```bash
python3 run_problem.py N1curl.ufl
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
