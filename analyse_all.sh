#!bin/bash


#for i in {2789..2817}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1873..1905}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1845..1864}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1865..1869}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1037..1045}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1870..1870}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1893..1905}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {1026..1035}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#xdone
#for i in {2832..2840}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done
#for i in {2731..2760}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done
#for i in {1799..1820}
#for i in 1822
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done
#for i in {2775..2781}
#do
#    ./drei_reader -a 2 $i | tee ${i}.log
#done
for i in 2764 2766 2767
do
    ./drei_reader -a 2 $i | tee ${i}.log
done



#do
#    sed -i 's|tsunamiA 12||' align/align_v2_${i}.dat 
#    sed -i 's|r152-scancal-tb21-1211.dat # missing|/mnt/pixeldata/a/r152-scancal-tb21-1211.dat|' align/align_v2_${i}.dat 
#    sed -i 's|r159-scancal-tb21-1210.dat #missing|/mnt/pixeldata/c/r159-scancal-tb21-1210.dat|' align/align_v2_${i}.dat 
#    sed -i 's|/home/cmspix/r4sclient/A/scm150-scancal1-tb21-2018-06-14-ia127-hold36.dat|/mnt/pixeldata/a/scm150-scancal1-tb21-2018-06-14-ia127-hold36.dat|' align/align_v2_${i}.dat 
#    ./drei_reader -a 2 $i | tee ${i}.log
#done

#for i in {2776..2781}
#do
#    cp align/align_v2_2775.dat align/align_v2_${i}.dat
#    git add align/align_v2_${i}.dat

#done
#git commit -m "added alignment 2775-2781, to be updated"    
#git push irene master
#for i in {2763..2768}
#do
 #  ./drei_reader -a 2 $i | tee ${i}.log

#done
    
