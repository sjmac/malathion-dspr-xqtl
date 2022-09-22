# malathion-dspr-xqtl
Code to analyze an X-QTL ("extreme QTL") mapping study of malathion resistance in a mixed population derived from DSPR (Drosophila Synthetic Population Resource) strains.

## Notes:

1. Shellscripts assume SLURM cluster scheduler.

2. Even if you use SLURM, header rows in the shellscripts will generally need to be changed to reflect your system.

3. Start with "get_bams_from_founders.txt" to get DSPR founder FASTQ data and align to reference yielding bam files. All based on using cluster.

4. Then move through "get_haplotype_frequencies.txt" to align X-QTL samples to genome (cluster), call SNPs (cluster), call haplotypes (cluster), compare haplotype frequencies between control and selected populations (local machine), and generate the figures from the paper (local machine).

## Walk through steps:

**File = get_bams_from_founders.txt**

  - Steps through getting the founder FASTQs from SRA, indexing reference, setting up cluster directory structure, and aligning reads to the reference. Assumes cluster for all steps.

**File = get_haplotype_frequencies.txt**

  - Steps through getting the X-QTL FASTQs from SRA, indexing reference, setting up cluster directory structure, and aligning reads to the reference for each sample. Then SNP/haplotype calling for both founders and X-QTL samples. All previous steps assumes using cluster. Subsequently, walks through executing malathion versus control haplotype tests, finding QTL peaks, and generating the plots in the paper, all assuming using a local machine. Note that if you want to simply recreate figures from the paper, can begin with "malathion_haps_genetic_1p5cM.txt.gz" (obtained from FigShare link associated with paper) and run through the code from "step 8" in the "get_haplotype_frequencies.txt" file.

## Brief explanation of files/folders:

**Folder = cluster_files_needed** (a series of files employed by shellscripts run on cluster)

  - **flymap.r6.txt** = Physical to genetic position conversion for Release 6 of the Dmel genome

**founders.txt** = Directs code to only examine the 8 DSPR pA founder strains

**readname.mapping.2acc.txt** = For founder assembly/variant calling for those founders with 2 SRA accessions

**readname.mapping.3acc.txt** = For founder assembly/variant calling for those founders with 3 SRA accessions

#

**Folder = cluster_r_scripts** (R scripts run via shellscripts on cluster)

**add_genetic_dist.R** = Add genetic distances to physical position files

**haplotyper_code_cM.R** = Actually runs the code to get haplotype frequencies for each pooled sample

**haplotyper_wrapper_cM.R** = A convenience wrapper to run the "haplotyper_code_cM.R" R script

**QC.freqs.R** = Executes some QC on the SNP calls (prior to haplotype calling)

#

**Folder = cluster_shellscripts** (shellscripts enabling analysis on the cluster)

**add_genetic_dist.sh** = Essentially just runs the "add_genetic_dist.R" R script

**bam2bcf.sh** = Runs bcftools to call SNPs

**founder_2acc_fq2bam.sh** = Read assembly (bwa, samtools) for founder strains with 2 SRA accessions

**founder_3acc_fq2bam.sh** = Read assembly (bwa, samtools) for founder strains with 3 SRA accessions

**fq2bam.sh** = Read assembly (bwa, samtools) for the pooled X-QTL samples

**haplotyper_cM.sh** = Essentially just runs the "haplotyper_wrapper_cM.R" R script

**QC.freqs.sh** = Largely just runs "QC.freqs.R" R script
