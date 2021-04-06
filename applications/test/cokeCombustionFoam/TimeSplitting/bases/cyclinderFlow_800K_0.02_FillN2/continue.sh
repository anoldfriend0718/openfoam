# decomposePar | tee -a run.log 
mpirun -np 8 $FOAM_USER_APPBIN/cokeCombustionFoam -parallel 2>&1 >> run.log 
reconstructPar  

# $FOAM_USER_APPBIN/dissolutionDBSFoam 2>&1 | tee -a run.log 