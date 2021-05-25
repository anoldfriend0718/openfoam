foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 


decomposePar | tee -a run.log 
mpirun -np 16 $FOAM_USER_APPBIN/icoDBSPermFoam -parallel 2>&1 >> run.log 
// reconstructPar | tee -a run.log 
