#!/bin/bash

read -p "GPU on or OFF, enter 0 for off and 1 for on?" GPU

RED='\033[0;31m'
NC='\033[0m'
echo " " 
echo " " 
echo "========================================================= " 
echo -e " ${RED}Compiling${NC} " 
rm bin/*
echo " " 
cd build/;cmake  ..;  make -j 8; cd ..
echo " " 
echo " " 
echo "========================================================= " 
echo " " 
echo -e "              ${RED}End of Compilation${NC}   " 
echo " " 
echo "========================================================= " 
echo " " 
echo -e "    removing ${RED} previous results ${NC} files ...     "
echo " " 
echo "--------------------------------------------------------- " 
echo " " 
echo " 1.removing Max.out " 
rm Max.out	
echo " " 
echo " 2.removing hd5 files" 
rm soln/*
echo " " 
echo " 3.removing meshSize0 " 
rm meshSize0
echo " " 
echo " " 
echo "========================================================= " 
echo " " 

echo -e "               ${RED}submitting job ....${NC}"
echo " " 
echo "--------------------------------------------------------- " 
echo " " 

read -p "enter the number of iteration for the run: " nn
printf "\n "
echo $nn

##for i in 1 4 9 16 25 36 49 64 81 100 121 144 169 196 225 256 289 324 361 
for((i=6;i<$nn;i=i+1)) 
do
echo "i = $i"
if [ $GPU -lt 1 ]; then
if [ $cluster == "bridges" ]; then
sbatch --ntasks $((i*i))  run_psc
elif [ $cluster == "stampede" ]; then
sbatch --ntasks $((i*i)) -N $(((i*i)/48+1))  run_tacc
elif [ $cluster == "crc" ]; then
sbatch --ntasks $((i*i)) -N $(((i*i)/28+1)) run_crc 
fi
else
sbatch --ntasks $((i*i))  run_gpu_psc.sh 
fi
done
echo " " 

sleep 1

#watch squeue -u shb105

echo "========================================================= " 
