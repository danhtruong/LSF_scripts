#!/bin/bash
#fasterq-dump script
#loadmodules
module load sratoolkit/2.9.4
module load parallel
#find directory
find $SCRATCH/opt/sra/*sra | parallel 'fasterq-dump -O $SCRATCH/opt/fastq {}'