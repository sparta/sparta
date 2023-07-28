#!/bin/bash -l

#PBS -N SPARTA_test
#PBS -m bea
#PBS -l select=1:ncpus=72:mpiprocs=72:mem=240GB
#PBS -l walltime=04:00:00
#PBS -q <queue>
#PBS -e sparta.err
#PBS -o sparta.out

module load hpcx/2.13.0
cd $PBS_O_WORKDIR

### Run program
mpirun -np 72 -bind-to socket -map-by socket ../../src/spa_kokkos_mpi_only < in.sparta

export SPARTA_PYTHON_TOOLS=/hpc/path/to/tools/pizza
python /path/to/pizza/log2txt.py log.sparta logdata.txt "Step" "Np" "v_mass_flow_x"
