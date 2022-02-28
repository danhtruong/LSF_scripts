#!/bin/bash/

# set up the software environment
module load samtools/1.14

echo "Processing file $f"

samtools sort -t CB -O BAM -o $f/outs/cellsorted_possorted_genome_bam.bam $f/outs/possorted_genome_bam.bam 



