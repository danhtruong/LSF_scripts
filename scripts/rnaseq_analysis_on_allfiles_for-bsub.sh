#! /bin/bash

for fq in $SCRATCH/opt/fastq/*.fastq
do
bsub -J rnaseq -W 3:00 -o output -e logs -q short -n 6 -M 32 -R rusage[mem=32] -N -B -u dtruong4@mdanderson.org
sh $HOME/scripts/rnaseq_analysis_on_input_file.sh $fq
sleep 1	# wait 1 second between each job submission
  
done

