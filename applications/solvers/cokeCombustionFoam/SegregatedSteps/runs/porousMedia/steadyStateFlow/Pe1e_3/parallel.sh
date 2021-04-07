foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 
# checkMesh  | tee  -a run.log 

# $FOAM_USER_APPBIN/icoDBSFoam  2>&1 | tee  -a run.log 

decomposePar | tee -a run.log 
mpirun -np 8 $FOAM_USER_APPBIN/icoDBSFoam -parallel 2>&1 | tee -a run.log 
reconstructPar | tee -a run.log 
