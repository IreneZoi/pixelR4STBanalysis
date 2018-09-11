#!bin/bash


for i in {1..1905}
do
file=align_$i.dat
file2=align_v2_$i.dat

    if [ -f $file ]; then
    cp  $file $file2
    echo "gainA " >> $file2
    echo "gainB " >> $file2
    echo "gainC " >> $file2
    echo "keA " >> $file2
    echo "keB " >> $file2
    echo "keC " >> $file2
    echo "beamEnergy " >> $file2
    echo "pitch " >> $file2
    echo "qL " >> $file2
    echo "qR " >> $file2
    echo "qLB " >> $file2
    echo "qRB " >> $file2
    echo "TsunamiA " >> $file2
    echo "TsunamiB " >> $file2
    echo "TsunamiC " >> $file2
    echo "dphcutA " >> $file2
    echo "dphcutB " >> $file2
    echo "dphcutC " >> $file2
    echo "dx3c " >> $file2
    fi
done
