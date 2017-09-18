#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="julia_mpi"
#SBATCH --output=julia_mpi.out
#SBATCH --mail-user=magnusnm@stud.ntnu.no
#SBATCH --mail-type=ALL
 
 
module load GCC
module load OpenMPI/2.0.2
make clean
make

echo
echo Running with 2 processors...
mpirun -np 2 ./julia_mpi -0.8 0.156
echo
echo Running with 4 processors...
mpirun -np 4 ./julia_mpi -0.8 0.156
echo
echo Running with 8 processors...
mpirun -np 8 ./julia_mpi -0.8 0.156
 
