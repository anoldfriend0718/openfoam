#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o reconstruct.log
#SBATCH -e err_reconstruct.log

. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

workerNum=60
caseDir=./
sampleRate=1

reconstruct $workerNum $caseDir $sampleRate

