#!/bin/bash/

# set up the software environment
module load python
module load samtools/1.14

echo "Processing file $f"

velocyto run10x -m $SCRATCH/hg38_rmsk.gtf -@ 6 --samtools-memory 64 $f $SCRATCH/refdata-gex-GRCh38-2020-A/genes/genes.gtf 

#velocyto run -b $f/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -o $f/velocyto -m $SCRATCH/hg38_rmsk.gtf $f/outs/cellsorted_possorted_genome_bam.bam $SCRATCH/refdata-gex-GRCh38-2020-A/genes/genes.gtf 

echo "Run complete"