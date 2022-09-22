#!/bin/bash
#SBATCH --job-name=callSNPs     # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sjmac@ku.edu     # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run (1 with >1 CPUS per task is multithreaded)
#SBATCH --cpus-per-task=1            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=2gb            # Memory request (default is 2Gb)
#SBATCH --time=48:00:00              # Time limit in hrs:min:sec (default for sjmac queue is 8 hrs)
#SBATCH --output=./log_files/callSNPs_%j.log    # Standard output and error logs
#SBATCH --array=1-5                  # Set number (N) of parallel jobs (1-N)

module load bwa/0.7.17
module load samtools/1.9
module load bcftools/1.9

ref="ref/dm6.fa"
bams="bams.txt"
dir1="process"
declare -a chrs=("chrX" "chr2L" "chr2R" "chr3L" "chr3R")
mychr=${chrs[$SLURM_ARRAY_TASK_ID - 1]}
bcftools mpileup -I -d 1000 -t $mychr -a "FORMAT/AD,FORMAT/DP" -f $ref -b $bams | bcftools call -mv -Ob > $dir1/calls.$mychr.bcf  
bcftools query -e'GT ="./."'  -e'QUAL<60' -f'%CHROM %POS %REF %ALT [ %AD{0} %AD{1}] [%GT]\n' $dir1/calls.$mychr.bcf | grep -v '\.' >$dir1/temp.$mychr.txt

