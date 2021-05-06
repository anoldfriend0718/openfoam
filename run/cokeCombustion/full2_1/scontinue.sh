#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o run1.log
#SBATCH -e err1.log

cpu=64
srun -n $cpu $FOAM_USER_APPBIN/cokeCombustionFoam2 -parallel
if [ $? -eq 0 ]; then
    echo "succeed to complete computations"
else
    echo "something wrong in computations"
fi

echo "start to reconstruct results"
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
