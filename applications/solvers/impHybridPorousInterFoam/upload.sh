#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

caseDir=./
runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/applications/solvers
foldername="$(basename `pwd`)"
syncFolderPath=$runspace/$foldername

echo "sync folder path on HPC: $syncFolderPath"

# papp_cloud rsync $caseDir/*.H          zc-m6:$syncFolderPath 
# papp_cloud rsync $caseDir/*.C          zc-m6:$syncFolderPath 
papp_cloud rsync $caseDir/*        bscc-a3:$syncFolderPath 
