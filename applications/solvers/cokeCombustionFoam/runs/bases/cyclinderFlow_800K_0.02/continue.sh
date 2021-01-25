# decomposePar | tee -a run.log 
mpirun -np 4 $FOAM_USER_APPBIN/cokeCombustionFoam -parallel 2>&1 >> run.log 
reconstructPar | tee -a run.log 

# $FOAM_USER_APPBIN/dissolutionDBSFoam 2>&1 | tee -a run.log 