#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 0-0:05:00
#SBATCH -J launch_psmc
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

SLIST=$1

for SAMPLE in $(cat ${SLIST})
do
 	sbatch vcf2fastq.sh ${SAMPLE}
done
