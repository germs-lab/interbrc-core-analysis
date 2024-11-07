#!/bin/bash --login
########## SBATCH Lines for Resource Request ##########

#SBATCH --time=40:00:00             # limit of wall clock time - how long the job will run (same as -t)
#SBATCH --nodes=26                   # number of different nodes - could be an exact number or a range of nodes (same as -N)
#SBATCH --ntasks=26                 # number of tasks - how many tasks (nodes) that you require (same as -n)
#SBATCH --cpus-per-task=1          # number of CPUs (or cores) per task (same as -c)
#SBATCH --mem-per-cpu=1G          # memory required per allocated CPU (or core)
#SBATCH --job-name 16S-QIIME-PROC       # you can give your job a name for easier identification (same as -J)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kristybr@msu.edu

########## Command Lines for Job Running ##########
module purge
module load Conda/3                                                 ### load necessary modules.
conda activate qiime2-amplicon-2024.10 

import_dir=/mnt/scratch/kristybr/raw_seqs                           # Directory for sequence imports (large number of small files)
dir=/mnt/research/EvansLab/Brandon/inter_BRC_microbiome/raw_seqs    # Directory for filtering and taxonomy 
cd $import_dir                                                      ### change to the directory where your code is located.

srun -n 26  job_name.sh                                             ### call your executable. (use srun instead of mpirun.)



scontrol show job $SLURM_JOB_ID                                     ### write job information to SLURM output file.
js -j $SLURM_JOB_ID                                                 ### write resource usage to SLURM output file (powertools command)



## Actual .sh script processing code 
## Perform genome alignments on subset of sequences from each BRC (CABBI, CBI, GLBRC, JBEI)(genome_alignment.sh)
#!/bin/bash
$bwa=/mnt/home/kristybr/bwa-0.7.17 ## Establish directory for bwa 
$ref=/mnt/home/kristybr/inter_BRC_merge/E_coli_genome_alignment ## Directory for reference genome and alignment outputs
## I randomly selected a .fastq file from each dataset to perform the alignment 


## GLBRC alignment 
cd /mnt/home/kristybr/inter_BRC_merge/GLBRC

$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna J10_S130_L001_R1_001.fastq > $ref/GLBRC_sample_R1.aln.sam 
$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna J10_S130_L001_R1_001.fastq > $ref/GLBRC_sample_R2.aln.sam 

## CABBI alignment 
cd /mnt/home/kristybr/inter_BRC_merge/CABBI

$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna 1001_C_2018_P3_N300_20180425_CABBI.R1.fastq > $ref/CABBI_sample_R1.aln.sam 
$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna 1001_C_2018_P3_N300_20180425_CABBI.R2.fastq > $ref/CABBI_sample_R2.aln.sam 

## JBEI alignment
cd /mnt/home/kristybr/inter_BRC_merge/JBEI 

$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna RHZA025_raw_1.fq > $ref/JBEI_sample_R1.aln.sam 
$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna RHZA025_raw_2.fq > $ref/JBEI_sample_R2.aln.sam 

## CBI alignment
cd /mnt/home/kristybr/inter_BRC_merge/CBI 

$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna 1016-Rhizo-16s-1-14-3_S48_L001_R1_001 > $ref/CBI_sample_R1.aln.sam 
$bwa/bwa mem $ref/GCF_000005845.2_ASM584v2_genomic.fna 1016-Rhizo-16s-1-14-3_S48_L001_R2_001 > $ref/CBI_sample_R2.aln.sam 

#Sam files were indexed and visualized with IGV genome view (figures attached) 
#GLBRC and CABBI have amplicons from 533BP to 785 BP.
#CBI and JBEI have amplicons from 514 BP to 805 BP.
# Use cutadapt to cut 19 BP from the forward and 20 BP from reverse reads for CBI/JBEI. All reads across all BRCs will then cover the same area. 

## Trim JBEI and CBI Reads with cutadapt (cutadapt.sh)
#!/bin/bash
module load Conda/3

conda activate cutadapt 

## Enter JBEI directory
cd /mnt/home/kristybr/inter_BRC_merge/JBEI

# Remove first 19BP from forward reads 
for file in *_R1.fq; do 
output_file="${file%.fq}_trimmed.fastq"
cutadapt -u 19 -o "$output_file" "$file"
done

# Remove last 20 BP from reverse reads.
for file in *_R2.fq; do 
output_file="${file%.fq}_trimmed.fastq"
cutadapt -u -20 -o "$output_file" "$file"
done

## Enter CBI directory 
cd /mnt/home/kristybr/inter_BRC_merge/CBI 

# Remove first 19BP from forward reads 
for file in *_R1_001.fastq; do 
output_file="${file%.fq}_trimmed.fastq"
cutadapt -u 19 -o "$output_file" "$file"
done

# Remove last 20 BP from reverse reads.
for file in *_R2_001.fastq; do 
output_file="${file%.fq}_trimmed.fastq"
cutadapt -u -20 -o "$output_file" "$file"
done


## Import FASTQ sequences into QIIME (import_seqs.sh)
#!/bin/bash 
# Define working directory                                         
module load Conda/3
conda activate qiime2-amplicon-2024.10 

import_dir=/mnt/scratch/kristybr/raw_seqs
cd $import_dir   

## Import demultiplexed, paired-end reads using PairedEndFastqManifestPhred33V2
# Convert UTF-8 .csv manifest file into a tab-delimited text file 
# sed 's/\,/\t/g' manifest_file.csv > manifest_file.txt

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path $import_dir/manifest_file.txt \
 --output-path interbrc_pe_demux.qza \
 --input-format PairedEndFastqManifestPhred33V2


## Filter and merge paired-end reads using dada2-denoise-paired (filter_seqs.sh)
#!/bin/bash
module load Conda/3
conda activate qiime2-amplicon-2024.10
dir=/mnt/research/EvansLab/Brandon/inter_BRC_microbiome/raw_Seqs
cd $dir 

qiime dada2 denoise-paired --i-demultiplexed-seqs interbrc_pe_demux.qza \ 
--p-trunc-len-f 250 \
--p-trunc-len-r 250 \
--p-n-threads 26 \
--o-representative-sequences interbrc_rep_seqs.qza \
--o-table interbrc_merged_feature_table.qza \
--o-denoising-stats dada2_denoising_stats.qza 

## Taxonomically classify the reads (classify_taxonomy.sh)
module load Conda/3
conda activate qiime2-amplicon-2024.10
dir=/mnt/research/EvansLab/Brandon/inter_BRC_microbiome/raw_Seqs
cd $dir 

qiime feature-classifier classify-sklearn --i-classifier classifier.qza \
--i-reads interbrc_rep_seqs.qza \
--o-classification interbrc_taxonomy.qza 

qiime metadata tabulate --m-input-file taxonomy.qza \
--o-visualization taxonomy.qzv 






