mpirun -np 8 $FOAM_USER_APPBIN/icoDBSPermFoam -parallel 2>&1 >> run.log 
reconstructPar | tee -a run.log 
