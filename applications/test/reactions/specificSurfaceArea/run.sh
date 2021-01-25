clear
wmake
$FOAM_USER_APPBIN/Test-SpecificSurfaceArea 2>&1 > run.log
foamToVTK -time 0
