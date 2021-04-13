clear
wmake
foamCleanTutorials
blockMesh
$FOAM_USER_APPBIN/Test-PMCCSpecificSurfaceArea 2>&1 > run.log
foamToVTK -time 0
