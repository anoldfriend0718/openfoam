#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o run.log
#SBATCH -e err.log

echo working directory: `pwd`

source   /public1/soft/openfoam/OpenFOAM7-fgl/OpenFOAM-7/etc/rebashrc
. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

blockMesh
# renumberMesh -overwrite -noFields  # not update the fields due to bug in the OpenFOAM on HPC
decomposePar

cpu=64
workerNum=$(($cpu-4)) #postProcessor 
caseDir=./
sampleRate=1 #postProcessing sampling rate
dataFolder=./postProcess
transverse_data_folder=./postProcess/transverseAveragedData/
imageFolder=./postProcess/images
animation_rate=20 #data interval 0.05 s

runWorkflow $cpu $workerNum $caseDir $sampleRate $dataFolder $transverse_data_folder $imageFolder $animation_rate


