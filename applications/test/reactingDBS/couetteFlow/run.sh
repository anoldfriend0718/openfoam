foamCleanTutorials 
logfile="run.log"
if [ -f $logfile ]; then
    rm $logfile
fi

touch run.log
blockMesh | tee -a $logfile
checkMesh | tee -a $logfile
$FOAM_USER_APPBIN/reactingDBSFoam 2>&1 > $logfile &

# $FOAM_USER_APPBIN/reactingDBSFoam 2>&1 | tee -a $logfile
# sleep 1 
# gnuplot ResidualPlots - 