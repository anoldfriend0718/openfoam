foamCleanTurtorial 
logfile="run.log"
if [ -f $logfile ]; then
    rm $logfile
fi

touch run.log
blockMesh | tee -a $logfile
checkMesh | tee -a $logfile
icoFoam 2>&1 | tee -a $logfile