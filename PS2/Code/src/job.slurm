#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="RPS_parallel"
#SBATCH --output=RPSparallel.out
#SBATCH --mail-user=magnusnm@stud.ntnu.no
#SBATCH --mail-type=ALL
 
 
module load GCC
module load OpenMPI/2.0.2
make clean
make

echo Starting run with 2 processors
mpirun -np 2 ./RPSparallel
 
echo Starting run with 4 processors
mpirun -np 4 ./RPSparallel

echo Starting run with 16 processors
mpirun -np 16 ./RPSparallel
