# decomposePar | tee -a run.log 
mpirun -np 8 $FOAM_USER_APPBIN/dbsReactiveThermalFoam2 -parallel 2>&1 >> run.log 
reconstructPar  
