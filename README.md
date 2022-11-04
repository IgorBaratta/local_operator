# Local Finite ELement Operator Benchmarks

Depends on:

- FFCx: The FEniCSx Form Compiler: [https://github.com/FEniCS/ffcx.git]
- For sum-factorization: [https://github.com/FEniCS/ffcx/tree/igor/tensor]

## Usage

```bash
python3 run.py --help
usage: run.py [-h] [--form_compiler {ffcx,ffc,tsfc}] [--scalar_type {double,float,_Float16,double _Complex,float _Complex}]
              [--problem {Laplacian,Mass,Elasticity,N1curl,Stokes}] [--conf CONF] [--degree DEGREE [DEGREE ...]] [--nrepeats NREPEATS]
              [--batch_size {None,1,2,4,8,16}] [--global_size GLOBAL_SIZE] [--action] [--mpi_size MPI_SIZE] [--cell_type {tetrahedron,hexahedron}]

Run local assembly benchmark.

optional arguments:
  -h, --help            show this help message and exit
  --form_compiler {ffcx,ffc,tsfc}
                        Form Compiler to use (default: ffcx)
  --scalar_type {double,float,_Float16,double _Complex,float _Complex}
                        Scalar type to use (default: double)
  --problem {Laplacian,Mass,Elasticity,N1curl,Stokes}
                        Problem to run (default: Laplacian)
  --conf CONF           Configuration file describing the compilers and flags. (default: compilers.yaml)
  --degree DEGREE [DEGREE ...]
                        Polynomial degree to evaluate the operators. (default: range(1, 4))
  --nrepeats NREPEATS   Number of times to run each experiment. (default: 3)
  --batch_size {None,1,2,4,8,16}
  --global_size GLOBAL_SIZE
                        Global number of dofs (assuming shared are dofs are duplicated). (default: 1000000.0)
  --action              Specify whether to run the problems with matrix free approach. (default: False)
  --mpi_size MPI_SIZE   The number of mpi processes to use. (default: 1)
  --cell_type {tetrahedron,hexahedron}
                        Cell type to use (default: tetrahedron
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

## Examples

### Matrix free Weighted Laplacian degrees 1-8

```bash
python3 run.py --problem Lagrange  --degree 1 2 3 4 5 6 7 8 --form_compiler=ffcx --action
```

### Single Precision Mass operator on 76 cores

```bash
python3 run.py --problem Mass --degree 1 2 3 4 5 6 7 8 --form_compiler=ffcx --action --mpi_size 76
```

## Output data description

Results:
output/{Problem}.txt

```bash
table = [machine,problem,compiler,version,flags,degree,form_compiler,scalar_type,batch_size,form_rank,cell_type,num_cells,time]
```
