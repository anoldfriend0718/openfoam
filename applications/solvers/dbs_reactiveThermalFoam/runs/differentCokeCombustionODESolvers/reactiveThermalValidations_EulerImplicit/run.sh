cd $WM_PROJECT_USER_DIR/applications/solvers/dbs_reactiveThermalFoam/
wmake
cd -
$FOAM_USER_APPBIN/dbsReactiveThermalFoam 2>&1 > run.log
