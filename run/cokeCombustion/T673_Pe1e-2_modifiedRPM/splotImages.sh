#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o plot.log
#SBATCH -e err_plot.log

. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

workerNum=60
dataFolder=./postProcess
imageFolder=./postProcess/images

plotImages $workerNum $dataFolder $imageFolder

