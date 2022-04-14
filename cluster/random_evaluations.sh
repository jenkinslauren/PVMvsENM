#!/bin/sh
#SBATCH -C [intel16|intel18]
#SBATCH -N 1 -c 1
#SBATCH -t 24:00:00
#SBATCH -o /mnt/research/TIMBER/PVMvsENM/QSTAT/random_evaluations.o
#SBATCH --mem 64G
#SBATCH -J random_evaluations

newgrp - TIMBER

module load GCC/8.3.0 OpenMPI/3.1.4 R/4.0.2

cd /mnt/research/TIMBER/PVMvsENM/code

Rscript random_evaluations.R

scontrol show job ${SLURM_JOB_ID}