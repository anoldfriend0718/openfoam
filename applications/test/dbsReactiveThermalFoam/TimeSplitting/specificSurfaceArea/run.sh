clear
wmake
foamCleanTutorials
blockMesh
$FOAM_USER_APPBIN/Test-ReactiveThermalCaseSpecificSurfaceArea 2>&1 > run.log
foamToVTK -time 0
