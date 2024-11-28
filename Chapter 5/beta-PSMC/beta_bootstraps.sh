#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH --array=1-20
#SBATCH -t 1-00:00:00
#SBATCH -J bPSMC
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools python3 python_ML_packages/3.9.5-gpu matplotlib/3.3.3-fosscuda-2020b

# run an analysis on the bootstrapped split files for each array
SAMPLE=$1

mkdir -p ${SAMPLE}_outfiles

/proj/sllstore2017021/nobackup/IVO/Beta-PSMC/beta-psmc -b -p 20*1 -N25 -t15 -r5 -B5 -o ${SAMPLE}_outfiles/b_${SLURM_ARRAY_TASK_ID}_${SAMPLE}.betapsmc.out ${SAMPLE}_infiles/bootstrapfiles/${SAMPLE}_split.psmcfa
python3 parse_psmc.py ${SAMPLE}_outfiles/b_${SLURM_ARRAY_TASK_ID}_${SAMPLE}.betapsmc.out beta_psmc bootstrap ${SLURM_ARRAY_TASK_ID} ${SAMPLE}
