#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o continue.log
#SBATCH -e err_continue.log

echo working directory: `pwd`

. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

cpu=64
workerNum=$(($cpu-4)) #postProcessor 
caseDir=./
sampleRate=10 #postProcessing sampling rate
dataFolder=./postProcess
imageFolder=./postProcess/images

runWorkflow $cpu $workerNum $caseDir $sampleRate $dataFolder $imageFolder








