#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o run.log
#SBATCH -e err.log
source   /public1/soft/openfoam/OpenFOAM7-fgl/OpenFOAM-7/etc/rebashrc
srun -n 64 $FOAM_USER_APPBIN/cokeCombustionFoam2 -parallel
