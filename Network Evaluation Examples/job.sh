#!/bin/bash

#SBATCH --time=24:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=32768M   # memory per CPU core
#SBATCH --mail-user=michael.bradshawiii@colorado.edu   # email address
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
which python3 > delete.txt
~/python3/bin/python3 michaels_test.py > output.txt 
