#!/bin/bash

runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/run/cokeCombustion
foldername="$(basename  "$(pwd)")"
syncFolderPath=$runspace/$foldername
echo "sync folder path: $syncFolderPath"

papp_cloud scp ../$foldername bscc-a3:$syncFolderPath 

