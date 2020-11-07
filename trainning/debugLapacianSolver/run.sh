#!/bin/zsh
blockMesh 2>&1 >run.log 
rm -rf postProcessing 
laplacianFoam 2>&1 >>run.log &
sleep 1
foamMonitor -l postProcessing/residuals/0/residuals.dat