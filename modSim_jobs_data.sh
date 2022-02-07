#!/bin/bash

#PBS -N modSim_5_March_final
#PBS -lnodes=1
#PBS -lwalltime=01:00:00

module load openmpi/gnu
module load R/3.2.3
module load eb
module load intel/2016b
module load fortran/intel

#export R_LIBS=$HOME/rpackages:$R_LIBS

cp -r "$HOME"/modSim "$TMPDIR"
cd "$TMPDIR"/modSim

Rscript --vanilla modSim_simulation_data.R iter

cp -r ./*.RDS "$HOME"/modSim/output



