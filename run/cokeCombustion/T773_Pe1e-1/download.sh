#!/bin/bash
cd ${0%/*} || exit 1    # Run from this directory
runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/run/cokeCombustion

foldername="$(basename  "$(pwd)")"
syncFolderPath=$runspace/$foldername
echo "sync folder path: $syncFolderPath"
papp_cloud rsync bscc-a3:$syncFolderPath/run*.log ./logs/
papp_cloud rsync bscc-a3:$syncFolderPath/err_*.log ./logs/
papp_cloud rsync bscc-a3:$syncFolderPath/postProcessing ./
papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/others ./postProcess/
papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/transverseAveragedData ./postProcess/
papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/animations ./postProcess/
papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/images/* ./postProcess/images/
papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/1.51.csv ./postProcess/
# papp_cloud rsync bscc-a3:$syncFolderPath/postProcess/*.csv ./postProcess/rawdata

# papp_cloud rsync bscc-a3:$syncFolderPath/0 ./
# papp_cloud rsync bscc-a3:$syncFolderPath/constant ./
# papp_cloud rsync bscc-a3:$syncFolderPath/system ./
# papp_cloud rsync bscc-a3:$syncFolderPath/s*.sh ./
# papp_cloud rsync bscc-a3:$syncFolderPath/*.py ./




