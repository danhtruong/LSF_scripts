#! /bin/bash
#adjust for naming 
#for SRA in $SCRATCH/opt/fastq/*R1_001.fastq

for SRA in $SCRATCH/opt/sra/*.sra
do
#base=$(awk -F'R1_001.fastq' '{print $1}' <<< $(basename $SRA))

fq1=$SCRATCH/opt/fastq/$(basename $SRA)_1.fastq
fq2=$SCRATCH/opt/fastq/$(basename $SRA)_2.fastq
bsub -J rnaseq -W 3:00 -o output -e logs -q short -n 28 -M 128 -R rusage[mem=128] -N -B -u dtruong4@mdanderson.org "bash $HOME/scripts/paired_rnaseq_analysis_on_input_file.sh $fq1 $fq2"
sleep 1	# wait 1 second between each job submission
  
done

