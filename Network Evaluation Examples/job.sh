#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=8   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32768M   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# arg 1 is the edge list for the network
# arg 2 is the disease gene sets
# arg 3 is the output dir

~/python3/bin/python3 evaluate_network.py $1 $2 $3 > $3/output.txt
