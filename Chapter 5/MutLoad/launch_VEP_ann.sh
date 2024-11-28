#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-00:05:00
#SBATCH -J VEP
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

# path to chromosome vcf file directory
SRCH=$1

# path to output directory
OUTDIR=$2

# txt with chromosome names
CHR=chr_list.txt

# parallelize by chromosome
for chr in $(cat ${CHR})
do
    sbatch VEP_ann.sh ${chr} ${SRCH} ${OUTDIR}
done


