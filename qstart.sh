#!/bin/bash

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

printf "\n "
echo $nn

sbatch  qrun_crc 
echo " " 

sleep 1

#watch squeue -u shb105

echo "========================================================= " 
