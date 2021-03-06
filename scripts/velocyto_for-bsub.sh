#! /bin/bash

for f in $SCRATCH/LPS_scRNA/velo/*; do
    if [ -d "$f" ]; then
        # $f is a directory
		echo "$f"
		export f=$f
		bsub < $HOME/LSF_scripts/velocyto_jobs.lsf
		sleep 1	# wait 1 second between each job submission
    fi
done

#samtools sort -t CB -O BAM -o LPS_scRNA/WDLS-4_analysis/outs/cellsorted_possorted_genome_bam.bam LPS_scRNA/WDLS-4_analysis/outs/possorted_genome_bam.bam