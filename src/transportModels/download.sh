#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory

caseDir=./
runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/src
foldername="$(basename `pwd`)"
syncFolderPath=$runspace/$foldername

echo "sync folder path on HPC: $syncFolderPath"

papp_cloud rsync bscc-a3:$syncFolderPath/* $caseDir      