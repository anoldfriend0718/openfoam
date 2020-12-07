foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 
checkMesh  | tee  -a run.log 

# decomposePar | tee -a $logfile
# mpirun -np 4 $FOAM_USER_APPBIN/reactingDBSFoam -parallel 2>&1 | tee -a $logfile
# reconstructPar | tee -a $logfile

$FOAM_USER_APPBIN/reactingDBSFoam 2>&1 | tee -a run.log 