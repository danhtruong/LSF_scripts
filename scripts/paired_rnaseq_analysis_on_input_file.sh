#!/bin/bash/

# This script takes a fastq file of RNA-Seq data, runs FastQC and outputs a counts file for it.
# USAGE: sh rnaseq_analysis_on_input_file.sh <name of fastq file>

# initialize a variable with an intuitive name to store the name of the input fastq file

fq1=$1
fq2=$2

# grab base of filename for naming outputs

base=$(basename $fq1 .subset.fastq)
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
align_out=$HOME/rnaseq/results/gsnap/${base}_Aligned.sortedByCoord.out.bam
counts=$HOME/rnaseq/results/counts/${base}_featurecounts.txt

# set up the software environment

module load fastqc
module load gmap/2020-06-04
module load samtools
module load subread

echo "Processing file $fq"

# Run FastQC and move output to the appropriate folder
fastqc $fq1 --outdir=$fastqc_out
fastqc $fq2 --outdir=$fastqc_out

# Run gsnap
gsnap -d $genome -D $genome_dir -t $cores --quality-protocol=sanger \
-M 2 -n 10 -B 2 -i 1 -N 1 -w 200000 -E 1 --pairmax-rna=200000 \
-A sam $fq1 $fq2 | samtools sort - | samtools view -bS - > $align_out

# Create BAM index
samtools index $align_out

# Count mapped reads
featureCounts -T $cores -s 2 -a $gtf -o $counts $align_out

