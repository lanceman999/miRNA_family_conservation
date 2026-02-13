#!/bin/bash

#SBATCH -J mirMachine                    # Job name
#SBATCH -A eande106                     # Allocation name
#SBATCH -p parallel                     # Partition/Queue name
#SBATCH -t 48:00:00                     # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                            # Number of nodes
#SBATCH -c 24                           # Number of cores


source activate mirmachine

# model proto for protosome - default is "combined" which was trained on miRNA from both proto and deutero - proto is supposed to be more accurate
# Mir-35 is not a family option for any model..... Mir-36 is present though 
base=$(basename $file)
strain=${base%%.*}

MirMachine.py --node Caenorhabditis --species $strain --genome $file --cpu 24 --model proto 
