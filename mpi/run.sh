#!/bin/bash
#SBATCH -N 2

DIR=$HOME/uoparallel-seismic-project/mpi

module load mpi
cd $DIR
./mt
