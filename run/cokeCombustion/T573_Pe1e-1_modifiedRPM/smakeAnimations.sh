#!/bin/bash
#SBATCH -p amd_256
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o animation.log
#SBATCH -e err_animation.log

. $WM_PROJECT_USER_DIR/utilities/scripts/postProcessFunctions

image_folder=./postProcess/images

cd $image_folder
animation_folder=../animations #relative to image_folder

rate=10

makeAnimation $animation_folder 'Temperature*.jpg'             4 temperature $rate &
makeAnimation $animation_folder 'porosity*.jpg'                4 porosity $rate &
makeAnimation $animation_folder 'O\$_2\$*.jpg'                 5 O2Conc $rate &
makeAnimation $animation_folder 'velocity-magnitude*.jpg '     5 velocityMagnitude $rate &
makeAnimation $animation_folder 'Reaction-Heat-Rate*.jpg '     6 heatRate $rate &

wait 

# ls Temperature*.jpg | sort -t '-' -k 4  -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > temperature.txt
 
# ffmpeg -y -r 1 -safe 0 -f concat  -i temperature.txt  -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/temperature.mp4
# ffmpeg  -y  -r 1 -safe 0 -f concat  -i temperature.txt -r 1   ../animations/temperature.gif

# rm temperature.txt

# ls porosity*.jpg | sort -t '-' -k 4 -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > porosity.txt
# ffmpeg -y -r 1 -safe 0 -f concat  -i porosity.txt  -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/porosity.mp4
# ffmpeg  -y  -r 1 -safe 0 -f concat  -i porosity.txt -r 1   ../animations/porosity.gif
# rm porosity.txt

# ls O\$_2\$*.jpg | sort -t '-' -k 5 -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > O2Conc.txt
# ffmpeg -y -r 1  -safe 0 -f concat  -i O2Conc.txt -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/O2Conc.mp4
# ffmpeg  -y  -r 1 -safe 0 -f concat  -i O2Conc.txt -r 1   ../animations/O2Conc.gif
# rm O2Conc.txt

# ls velocity-magnitude*.jpg | sort -t '-' -k 5 -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > velocitymagnitude.txt
# ffmpeg -y -r 1 -safe 0 -f concat  -i velocitymagnitude.txt  -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/velocitymagnitude.mp4
# ffmpeg  -y  -r 1 -safe 0 -f concat  -i velocitymagnitude.txt -r 1   ../animations/velocitymagnitude.gif
# rm velocitymagnitude.txt


# ls Reaction-Heat-Rate*.jpg | sort -t '-' -k 6 -n | sed 's:\ :\\\ :g'| sed 's/^/file /'  > ReactionHeatRate.txt
# ffmpeg -y -r 1 -safe 0 -f concat  -i velocitymagnitude.txt  -r 5 -vcodec libx264  -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2"  ../animations/ReactionHeatRate.mp4
# ffmpeg  -y  -r 1 -safe 0 -f concat  -i ReactionHeatRate.txt -r 1   ../animations/ReactionHeatRate.gif
# rm ReactionHeatRate.txt

