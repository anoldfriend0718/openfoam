# decomposePar | tee -a run.log 
mpirun -np 8 $FOAM_USER_APPBIN/icoDBSFoam  -parallel 2>&1 >> run.log 
reconstructPar  
