#!/bin/bash
######################################################
# qsub arguments
######################################################
#SBATCH --job-name=tl-%seed%
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=long
#SBATCH --mem-per-cpu=3000MB
######################################################
# lib 
######################################################
######################################################
# Collect arguments
######################################################
source /home/yinc/miniconda3/bin/activate p3.8


for x in $(seq 1 10);do
    python3 test-%seed%-${x}.py > /data/yinc/ch2oo-project/MD-simulations/13877/output%seed%-${x}.out
    python3 check_last_frame.py -i nve_seed%seed%-${x}.dcd
    return_code=$?
    if [ $return_code -eq 0 ];then
       break
    fi
done


