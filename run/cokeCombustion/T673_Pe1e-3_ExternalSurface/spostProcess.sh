#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o postProcess.log
#SBATCH -e err_postProcess.log


. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

data_folder=./postProcess
transverse_data_folder=./postProcess/transverseAveragedData/
cpu=64

batchPostProcess $data_folder $transverse_data_folder $cpu