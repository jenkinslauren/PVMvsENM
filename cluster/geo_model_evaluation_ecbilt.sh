#!/bin/sh
#SBATCH -C [intel16|intel18]
#SBATCH -N 1 -c 1
#SBATCH -t 24:00:00
#SBATCH -o /mnt/research/TIMBER/PVMvsENM/QSTAT/geo_model_evaluation_ecbilt.o
#SBATCH --mem 64G
#SBATCH -J geo_model_evaluation_ecbilt

newgrp - TIMBER

module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

cd /mnt/research/TIMBER/PVMvsENM/code

Rscript geo_model_evaluation.R 'ecbilt'

scontrol show job ${SLURM_JOB_ID}