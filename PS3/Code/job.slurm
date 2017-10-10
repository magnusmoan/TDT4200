#!/bin/sh
#SBATCH --partition=TRAINING
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=12000
# Memory per node specification is in MB. It is optional.
# The default limit is 3000MB per core.
#SBATCH --job-name="RPS_pthread"
#SBATCH --output=foo.out
#SBATCH --mail-user=magnusnm@stud.ntnu.no
#SBATCH --mail-type=ALL

module load GCC
make clean
make pth

./pthreadRPS.out 2
