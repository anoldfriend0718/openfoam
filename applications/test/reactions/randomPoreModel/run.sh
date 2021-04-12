clear
wmake
$FOAM_USER_APPBIN/Test-RandomPoreModel 2>&1 > run.log
foamToVTK -time 0
