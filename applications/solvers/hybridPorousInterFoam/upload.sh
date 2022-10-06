#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

# caseDir=./
# runspace=/public1/home/scg1301/OpenFOAM/scg1301-7/applications/solvers
# # runspace=/public4/home/sc56333/OpenFOAM/sc56333-7/applications/solvers
# foldername="$(basename `pwd`)"
# syncFolderPath=$runspace/$foldername

# echo "sync folder path on HPC: $syncFolderPath"

# # papp_cloud rsync $caseDir/*          zc-m6:$syncFolderPath 
# # papp_cloud rsync $caseDir/0_org          zc-m6:$syncFolderPath 
# # papp_cloud rsync $caseDir/constant   zc-m6:$syncFolderPath
# # papp_cloud rsync $caseDir/system     zc-m6:$syncFolderPath
# # papp_cloud rsync $caseDir/s*.sh      zc-m6:$syncFolderPath
# # papp_cloud rsync $caseDir/*.py       bscc-a:$syncFolderPath

# papp_cloud rsync $caseDir/*.H          bscc-a:$syncFolderPath 
# papp_cloud rsync $caseDir/*.C          bscc-a:$syncFolderPath 


caseDir=./
runspace=/public1/home/scg1301/OpenFOAM/scg1301-7/applications/solvers
foldername="$(basename `pwd`)"
syncFolderPath=$runspace/$foldername

echo "sync folder path on HPC: $syncFolderPath"

# papp_cloud rsync $caseDir/*.H          zc-m6:$syncFolderPath 
# papp_cloud rsync $caseDir/*.C          zc-m6:$syncFolderPath 
papp_cloud rsync $caseDir/*        zc-m6:$syncFolderPath 
