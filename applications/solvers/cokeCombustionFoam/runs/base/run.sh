#
cd  ../../
wmake
cd -
$FOAM_USER_APPBIN/cokeCombustionFoam 2>&1 > run.log
