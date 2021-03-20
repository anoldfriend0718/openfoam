cd $WM_PROJECT_USER_DIR/applications/solvers/dbs_conjugateHeatTransfer/
wmake
cd -
$FOAM_USER_APPBIN/dbsConjugateHeatTransfer 2>&1 > run.log
