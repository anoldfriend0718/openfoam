#!/bin/bash

runspace=/public1/home/sc90898/OpenFOAMWorkspace/sc90898-7/run/cokeCombustion
syncfile=residuals*.jpg
foldername="$(basename  "$(pwd)")"
syncFilePath=$runspace/$foldername/$syncfile
echo "sync file path: $syncFilePath"

while true

do

papp_cloud rsync bscc-a3:$syncFilePath ./

sleep 10

done

