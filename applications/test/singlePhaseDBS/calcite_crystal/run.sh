foamCleanTutorials | tee run.log 
blockMesh  | tee -a run.log 
checkMesh  | tee  -a run.log 
rm 0 -rf  
cp 0.orig 0 -rf 
$FOAM_USER_APPBIN/dbs_pimpleFoamCi 2>&1 | tee  -a run.log 