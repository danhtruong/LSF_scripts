#!/bin/bash/

# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file

fq=$1

# grab base of filename for naming outputs

base=$(basename $fq .subset.fastq)
echo "Sample name is $base"           

# specify the number of cores to use

cores=28

# directory with genome reference FASTA and index files + name of the gene annotation file

genome=hg38
genome_dir=$SCRATCH/hg38/
gtf=$SCRATCH/hg38/GCF_000001405.39_GRCh38.p13_genomic.gtf

# make all of the output directories
# The -p option means mkdir will create the whole path if it 
# does not exist and refrain from complaining if it does exist

mkdir -p $HOME/rnaseq/results/fastqc/
mkdir -p $HOME/rnaseq/results/gsnap/
mkdir -p $HOME/rnaseq/results/counts/

# set up output filenames and locations

fastqc_out=$HOME/rnaseq/results/fastqc/
align_out_1=$HOME/rnaseq/results/gsnap/SRR7014594.sra_1.fastq_Aligned.sortedByCoord.out.bam
align_out_2=$HOME/rnaseq/results/gsnap/SRR7014594.sra_2.fastq_Aligned.sortedByCoord.out.bam
counts=$HOME/rnaseq/results/counts/SRR7014594.sra.fastq_featurecounts.txt

# set up the software environment

module load subread

# Count mapped reads
featureCounts -T $cores -s 2 -a $gtf -o $counts $align_out_1 $align_out_2