#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

caseDir=./
runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/run/digitalRockPerm
foldername="$(basename  "$(pwd)")"
syncFolderPath=$runspace/$foldername

echo "sync folder path on HPC: $syncFolderPath"

papp_cloud rsync $caseDir/0          bscc-a3:$syncFolderPath 
papp_cloud rsync $caseDir/constant   bscc-a3:$syncFolderPath
papp_cloud rsync $caseDir/system     bscc-a3:$syncFolderPath
papp_cloud rsync $caseDir/s*.sh      bscc-a3:$syncFolderPath
papp_cloud rsync $caseDir/*.py       bscc-a3:$syncFolderPath



