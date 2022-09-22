#!/bin/bash
#SBATCH --job-name=founder_bwa_2acc  # Name of job
#SBATCH --mail-type=BEGIN,END,FAIL   # Set when I get emails (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sjmac@ku.edu     # Email address
#SBATCH --partition=sjmac            # Either sixhour or sjmac
#SBATCH --nodes=1                    # Number of nodes to run on
#SBATCH --ntasks=1                   # Number of tasks to run (1 with >1 CPUS per task is multithreaded)
#SBATCH --cpus-per-task=4            # Number of CPUs per task (>1 = multithreaded)
#SBATCH --mem-per-cpu=2gb            # Memory request (default is 2Gb)
#SBATCH --time=6:00:00               # Time limit in hrs:min:sec (default for sjmac queue is 8 hrs)
#SBATCH --output=./run_log_files/found_bwa_2acc_%j.log     # Standard output and error log
#SBATCH --array=1-3                  # Set number (N) of parallel jobs (1-N)

module load bwa/0.7.17
module load samtools/1.9
module load picard/2.20.3

# ASSUMES specific directory structure
ref="ref/dm6.fa"
dir1="process"
dir2="bam"
files="readname.mapping.2acc.txt"

# ASSUMES "readname.mapping.2acc.txt" has format
# Col1   = Desired name of sample
# Col2/3 = Read1/Read2 of sequencing replicate "a"
# Col4/5 = Read1/Read2 of sequencing replicate "b"
shortname=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f1`
F1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f2`
R1=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f3`
F2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f4`
R2=`head -n $SLURM_ARRAY_TASK_ID $files | tail -n 1 | cut -f5`

bwa mem -t 4 -M $ref $F1 $R1 | samtools view -bS - > $dir1/$shortname.1.temp.bam
bwa mem -t 4 -M $ref $F2 $R2 | samtools view -bS - > $dir1/$shortname.2.temp.bam
samtools merge -f $dir1/$shortname.temp.bam $dir1/$shortname.1.temp.bam $dir1/$shortname.2.temp.bam
samtools sort $dir1/$shortname.temp.bam -o $dir2/$shortname.bam
picard AddOrReplaceReadGroups I=$dir2/$shortname.bam O=$dir2/$shortname.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$shortname RGSM=$shortname VALIDATION_STRINGENCY=LENIENT
samtools index $dir2/$shortname.RG.bam

