foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 
checkMesh  | tee  -a run.log 

decomposePar | tee -a run.log 
mpirun -np 8 $FOAM_USER_APPBIN/dbsConjugateHeatTransfer -parallel 2>&1 >> run.log 
reconstructPar #| tee -a run.log 

# $FOAM_USER_APPBIN/dissolutionDBSFoam 2>&1 | tee -a run.log 