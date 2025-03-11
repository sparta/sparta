#!/bin/bash
#SBATCH --time=01:00:00 
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=couple_test
#SBATCH --partition=low
#SBATCH --output=std.out
#SBATCH --error=std.err

mpirun spa.exe -in in.couple
