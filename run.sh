#!/bin/bash -l
#SBATCH --job-name=mpi_job_test      # Job name
#SBATCH --cpus-per-task=1            # Number of cores per MPI task 
#SBATCH --nodes=1                    # Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=8          # Maximum number of tasks on each node
#SBATCH --output=mpi_test_%j.log     # Path to the standard output and error files relative to the working directory
#SBATCH -p small

spack env activate ffcx
spack load cmake
python3 run.py --problem Elasticity  --degree 2 --form_compiler=ffcx --action --global_size 10000000 --output_file=output/elasticity_2_action.csv
python3 run.py --problem Elasticity  --degree 6 --form_compiler=ffcx --action --global_size 10000000 --output_file=output/elasticity_6_action.csv
