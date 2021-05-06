#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o reconstruct.log
#SBATCH -e reconstruct.log

cpu=64
source $WM_PROJECT_USER_DIR/.venv/bin/activate
caseDir=./
pyscipt=$WM_PROJECT_USER_DIR/utilities/postProcess/pyResconstruct.py
fieldNames='["U","T","p","rho","O2","CO2","eps","coke","cokeRectionRate","Qdot"]'
timeNames=all
sampleRate=10
workerNum=$cpu
overWrite=false
dataFolder=postProcess

if $overWrite; then
    python3  $pyscipt -c $caseDir -t $timeNames -f $fieldNames -n $workerNum -r $sampleRate -s $dataFolder -w
else
    python3  $pyscipt -c $caseDir -t $timeNames -f $fieldNames -n $workerNum -r $sampleRate -s $dataFolder
fi
