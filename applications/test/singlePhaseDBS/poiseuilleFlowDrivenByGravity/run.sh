foamCleanTutorials 
logfile="run.log"
if [ -f $logfile ]; then
    rm $logfile
fi

touch run.log
blockMesh | tee -a $logfile
checkMesh | tee -a $logfile
$FOAM_USER_APPBIN/icoGravityFoam 2>&1 | tee -a $logfile