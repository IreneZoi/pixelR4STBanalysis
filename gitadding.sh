#!bin/bash

for i in {1026..1035}
do
#    cp /mnt/pixeldata/a/align_${i}.dat  align/align_v2_${i}.dat
    git add align/align_v2_${i}.dat

done
git commit -m "added alignment 1026-1035, up to date"    
git push irene master
    
