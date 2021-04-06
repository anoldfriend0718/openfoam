cd $WM_PROJECT_USER_DIR/applications/solvers/cokeCombustionFoam/TimeSplitting 
wmake
cd -
$FOAM_USER_APPBIN/cokeCombustionFoam 2>&1 > run.log
