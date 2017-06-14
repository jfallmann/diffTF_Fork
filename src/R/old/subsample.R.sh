#! /bin/bash
#BSUB -J subsample
#BSUB -M 100000
#BSUB -n 125
#BSUB -q medium_priority
#BSUB -o /scratch/carnold/CLL/TF_act_downsampling0.125/output/Logs_and_Benchmarks/subsample_output.txt
#BSUB -e /scratch/carnold/CLL/TF_act_downsampling0.125/output/Logs_and_Benchmarks/subsample_error.txt

echo "I, ${USER}, am running on host ${HOSTNAME}"

echo "STARTING SNAKEMAKE MASTER PROCESS"

Rscript /g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/subsample.R

echo "FINISHED SNAKEMAKE MASTER PROCESS"

# Submit with bsub < BSUB_TEMPLATE.sh
# Check status with bjobs
# Check output during run with bpeek JOBID
