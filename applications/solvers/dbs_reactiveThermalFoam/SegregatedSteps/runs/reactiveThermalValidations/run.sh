cd $WM_PROJECT_USER_DIR/applications/solvers/dbs_reactiveThermalFoam/SegregatedSteps
wmake
cd -
$FOAM_USER_APPBIN/dbsReactiveThermalFoam2 2>&1 > run.log