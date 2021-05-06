#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o reconstruct.log
#SBATCH -e reconstruct.log


data_folder=./postProcess
image_folder=./postProcess/images
animation_folder=./postProcess/animations
fields='["eps","UNorm","O2Conc","Qdot","T"]'
times='all'
script=${WM_PROJECT_USER_DIR}/utilities/postProcess/pyFigure.py
workerNum=64

source $WM_PROJECT_USER_DIR/.venv/bin/activate
echo "start to plot images"
python3 $script -f $fields -d $data_folder -s $image_folder -t $times -n $workerNum

if [ $? -eq 0 ]; then
    echo "succeed to complete plotting "
    echo "start to make animations"
    if [ -f $image_folder ]; then
        echo "File \"/path/to/file\" exists"
    fi

    if [ -e $image_folder ]; then
        echo "Path $image_folder exists"
    else
        mkdir $image_folder
    fi
    
    
    cd $image_folder
    ffmpeg -y  -framerate 1 -pattern_type glob -i "porosity*.jpg"   -r 10  ../animations/porosity.gif
    ffmpeg -y  -framerate 1 -pattern_type glob -i "velocitymagnitude*.jpg"  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 10  ../animations/velocitymagnitude.gif
    ffmpeg -y  -framerate 1 -pattern_type glob -i "Temperature*.jpg"   -r 10  ../animations/temperature.gif
    ffmpeg -y  -framerate 1 -pattern_type glob -i "ReactionHeatRate*.jpg"   -r 10  ../animations/reactionHeat.gif
    ffmpeg -y  -framerate 1 -pattern_type glob -i "O\$_2\$conc*.jpg"   -r 10  ../animations/O2conc.gif
else
    echo "something wrong in plotting"
fi