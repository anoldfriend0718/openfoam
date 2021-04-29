#!/bin/bash
source $WM_PROJECT_USER_DIR/.venv/bin/activate
caseDir=$WM_PROJECT_USER_DIR/applications/solvers/cokeCombustionFoam/SegregatedSteps/runs/complicatedPorousMedia/combustions/tiny/convectionLimited
pyscipt=$WM_PROJECT_USER_DIR/utilities/postProcess/pyResconstruct.py
fieldNames='["U","T","p","rho","O2","CO2","eps","coke","cokeRectionRate","Qdot"]'
timeNames=['0.1']
sampleRate=1
workerNum=30
overWrite=false
dataFolder=postProcess

if $overWrite; then
    python3  $pyscipt -c $caseDir -t $timeNames -f $fieldNames -n $workerNum -r $sampleRate -s $dataFolder -w
else
    python3  $pyscipt -c $caseDir -t $timeNames -f $fieldNames -n $workerNum -r $sampleRate -s $dataFolder 
fi
