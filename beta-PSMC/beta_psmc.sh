#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1-00:00:00
#SBATCH -J bPSMC
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error


module load bioinfo-tools python3 python_ML_packages/3.9.5-gpu matplotlib/3.3.3-fosscuda-2020b

SAMPLE=$1

## first concatenate all the fasta files
for chrom in $(cat chrm_auto_intervals.txt | cut -f 1)
do
 	zcat ${SAMPLE}_infiles/${chrom}_${SAMPLE}.fa.gz >> ${SAMPLE}_infiles/${SAMPLE}.fa
done

gzip ${SAMPLE}_infiles/${SAMPLE}.fa

rm -rf ${SAMPLE}_infiles/NC_*.fa.gz

# prep the psmc input
/proj/sllstore2017021/nobackup/IVO/Beta-PSMC/utils/fq2psmcfa -q20 ${SAMPLE}_infiles/${SAMPLE}.fa.gz > ${SAMPLE}_infiles/${SAMPLE}.psmcfa

# run the psmcfa splitter 
## splitfa is a script that comes with the beta-psmc package that splits the input sequence sinto shorter segments for random bootstrapping
mkdir -p ${SAMPLE}_infiles/bootstrapfiles
/proj/sllstore2017021/nobackup/IVO/Beta-PSMC/utils/splitfa ${SAMPLE}_infiles/${SAMPLE}.psmcfa > ${SAMPLE}_infiles/bootstrapfiles/${SAMPLE}_split.psmcfa
mkdir -p ${SAMPLE}_outfiles

# submit an array job for bootstrapping
sbatch --array=1-20 beta_bootstraps.sh ${SAMPLE}

# and submit the regular psmc too # this is optional
#sbatch psmc.sh ${SAMPLE}

# and then run the analysis
/proj/sllstore2017021/nobackup/IVO/Beta-PSMC/beta-psmc -p 20*1 -N25 -t15 -r5 -B5 -o ${SAMPLE}_outfiles/${SAMPLE}.betapsmc_main.out ${SAMPLE}_infiles/${SAMPLE}.psmcfa
# parse the output
python3 parse_psmc.py ${SAMPLE}_outfiles/${SAMPLE}.betapsmc_main.out beta_psmc main 0 ${SAMPLE}