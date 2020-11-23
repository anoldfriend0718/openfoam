foamCleanTutorials

logfile="run.log"
if [ -f $logfile ]; then
    rm $logfile
fi

touch run.log
blockMesh | tee -a $logfile
checkMesh | tee -a $logfile

decomposePar | tee -a $logfile
mpirun -np 4 icoFoam -parallel 2>&1 | tee -a $logfile
reconstructPar | tee -a $logfile

