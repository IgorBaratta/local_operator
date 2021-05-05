# Benchmark FFCx

Depends on:
- ffcx -  branch https://github.com/FEniCS/ffcx/tree/igor/loop-invariant


## How to run:
### Lagrange - Stiffness Matrix - problem 0
```bash
python3 run_problem.py 0
```
### N1curl - curl curl - problem 1
```bash
python3 run_problem.py 1
```

### DG - ?? - problem 2
```bash
python3 run_problem.py 2
```

Results:
{Name}+degre.txt


## Add more compilers or flags:
Append compiler names and flags to lines 10 and 11 of [run_problem](https://github.com/IgorBaratta/benchmark_function/blob/main/run_problem.py)


## Plotting data:
```bash
python3 graph.py degree
python3 graph.py 1
python3 graph.py 2
```

## Data description
```
table = [method, compiler name, compiler flags, polyonomial degree, number of cells, local assemble time]
```
