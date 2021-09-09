#!/bin/bash
#SBATCH --time=00:00:30               # Run time (days-hh:mm:ss) - (max 7days)
#SBATCH --partition=batch             # Submit to queue/partition named batch
#SBATCH --nodes=5
#SBATCH --job-name=myTest
#SBATCH --cpus-per-task=10

module load gcc openmpi
module load gcc
make final