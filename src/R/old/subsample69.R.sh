#! /bin/bash
#BSUB -J subsample
#BSUB -M 100000 
#BSUB -n 40
#BSUB -q medium_priority
#BSUB -o /scratch/carnold/CLL/TF_act/output/Logs_and_Benchmarks/subsample_Ivan_89output.txt
#BSUB -e /scratch/carnold/CLL/TF_act/output/Logs_and_Benchmarks/subsample_Ivan_89error.txt

echo "I, ${USER}, am running on host ${HOSTNAME}"

echo "STARTING SNAKEMAKE MASTER PROCESS"

Rscript /g/scb2/zaugg/carnold/Projects/PWM_Ivan/src/R/subsample.IVan69.R

echo "FINISHED SNAKEMAKE MASTER PROCESS"

# Submit with bsub < BSUB_TEMPLATE.sh
# Check status with bjobs
# Check output during run with bpeek JOBID
