foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 
# checkMesh  | tee  -a run.log 



decomposePar | tee -a run.log 

mpirun -np 8 $FOAM_USER_APPBIN/cokeCombustionFoam2 -parallel 2>&1 >> run.log 
./reconstruct.sh 

