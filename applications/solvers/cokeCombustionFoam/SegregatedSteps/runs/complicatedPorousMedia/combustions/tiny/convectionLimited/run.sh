cd $WM_PROJECT_USER_DIR/applications/solvers/cokeCombustionFoam/SegregatedSteps
wmake
cd -
$FOAM_USER_APPBIN/cokeCombustionFoam2 2>&1 > run.log