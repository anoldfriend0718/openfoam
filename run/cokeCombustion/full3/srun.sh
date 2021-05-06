#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o run.log
#SBATCH -e err.log
source   /public1/soft/openfoam/OpenFOAM7-fgl/OpenFOAM-7/etc/rebashrc
cpu=64
blockMesh
renumberMesh -overwrite -noFields  # not update the fields due to bug in the OpenFOAM on HPC
decomposePar
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

