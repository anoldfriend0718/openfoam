#!/bin/bash
# runspace=/public4/home/sc56333/OpenFOAM/sc56333-7/applications/solvers

# foldername="$(basename `pwd`)"
# folderpath=$runspace/$foldername

# echo HPC directory: $folderpath
# papp_cloud ssh bscc-a



runspace=/public1/home/scg1301/OpenFOAM/scg1301-7/applications/solvers
foldername="$(basename `pwd`)"
folderpath=$runspace/$foldername

echo HPC directory: $folderpath
papp_cloud ssh zc-m6