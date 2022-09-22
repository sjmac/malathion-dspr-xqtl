#!/bin/bash
#SBATCH --job-name=add_cM_info       # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sjmac@ku.edu     # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run (1 with >1 CPUS per task is multithreaded)
#SBATCH --cpus-per-task=4            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=2gb            # Memory request (default is 2Gb)
#SBATCH --time=6:00:00               # Time limit in hrs:min:sec (default for sjmac queue is 8 hrs)
#SBATCH --output=./log_files/add_cM_%j.log    # Standard output and error logs

module load R/4.0
# TELL R where to find libraries
# because 'add_genetic_dist.R' called below calls 'library(tidyverse)'
export R_LIBS_USER=/panfs/pfs.local/work/sjmac/sjmac/R

Rscript scripts/add_genetic_dist.R "FREQ_SNPs.txt"

