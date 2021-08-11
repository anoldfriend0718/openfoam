#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o run1.log
#SBATCH -e err1.log

echo working directory: `pwd`

. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

cpu=64
workerNum=$(($cpu-4)) #postProcessor 
caseDir=./
sampleRate=1 #postProcessing sampling rate
dataFolder=./postProcess
transverse_data_folder=./postProcess/transverseAveragedData/
imageFolder=./postProcess/images
animation_rate=20 #data interval 0.05 s


runWorkflow $cpu $workerNum $caseDir $sampleRate $dataFolder $transverse_data_folder $imageFolder $animation_rate











