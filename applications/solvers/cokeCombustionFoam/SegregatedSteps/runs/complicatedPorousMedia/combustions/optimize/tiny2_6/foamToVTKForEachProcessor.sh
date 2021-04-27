#!/bin/bash                        
## 

dir=~/OpenFOAM/anoldfriend-7/applications/solvers/cokeCombustionFoam/SegregatedSteps/runs/complicatedPorousMedia/combustions/optimize/tiny2_6
for ((i=0; i<=7; i++))
do
   cd $dir/processor$i
   foamToVTK
done
