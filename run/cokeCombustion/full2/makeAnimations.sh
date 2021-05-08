#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o animation.log
#SBATCH -e err_animation.log

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

    if [ -e $animation_folder ]; then
        echo "Animation path exists: $animation_folder"
    else
        echo "Create animation path: $animation_folder"
        mkdir $animation_folder
    fi

    # cd $image_folder
    
    # ffmpeg -y  -framerate 1 -pattern_type glob -i "porosity*.jpg"   -r 10  ../animations/porosity.gif
    # ffmpeg -y  -framerate 1 -pattern_type glob -i "velocitymagnitude*.jpg"  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -r 10  ../animations/velocitymagnitude.gif
    # ffmpeg -y  -framerate 1 -pattern_type glob -i "Temperature*.jpg"   -r 10  ../animations/temperature.gif
    # ffmpeg -y  -framerate 1 -pattern_type glob -i "ReactionHeatRate*.jpg"   -r 10  ../animations/reactionHeat.gif
    # ffmpeg -y  -framerate 1 -pattern_type glob -i "O\$_2\$conc*.jpg"   -r 10  ../animations/O2conc.gif

    # echo "succeed to make animations"

else
    echo "something wrong in plotting"
fi

# cd $image_folder

# ls O\$_2\$-conc*.jpg | sort -t '-' -k 4 -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > O2Conc.txt
 
# ffmpeg -y -r 1 -f concat -safe 0  -i O2Conc.txt  -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/O2Conc.mp4
# ffmpeg  -y  -r 1 -f concat -safe 0  -i O2Conc.txt -r 1   ../animations/O2Conc.gif

# rm O2Conc.txt
# ffmpeg  -y  -framerate 1  -i concat:$CONCAT    -r 10   ../animations/temperature.gif
