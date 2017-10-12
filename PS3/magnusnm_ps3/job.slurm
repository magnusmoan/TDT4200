#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12000
#SBATCH --job-name="RPS_pthread"
#SBATCH --output=threading_%.out
#SBATCH --mail-user=magnusnm@stud.ntnu.no
#SBATCH --mail-type=ALL

module load GCC
make clean

make omp
echo Starting omp
./ompRPS.out 16
echo Done with omp

make pth
echo Starting pthread
./pthreadRPS.out 16
echo Done with pthread

