#!/bin/bash
#SBATCH -J cbtest       # Job Name
#SBATCH -o DNA.out%j    # Output and error file name (%j expands to jobID)
#SBATCH -n 1           # Total number of mpi tasks requested
#SBATCH -p serial  # Queue (partition) name -- normal, development, etc.
#SBATCH -t 02:00:00     # Run time (hh:mm:ss) - 1.5 hours
#SBATCH -A GKIMP     # Project name
###ibrun ../bin/dna           # Run the MPI executable 
module load python
python cb_test_script.py

