#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-01:00:00
#SBATCH -J out_align
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools minimap2 samtools bcftools

# path to chromosome vcf file
chr=$1

# path to population level vcf files
SRCH=$2

# path to outgroup vcf
OUTG=$3

# path to output directory
OUTDIR=$4

# output name
SMP=$4

# get the chromosome number to restrict merging the vcf files
## this will depend heavily on how your files are named; don't forget to tweak this part
## the end result should be something like NC_041754.1
chr_pos=$(echo ${chr} | cut -d "." -f 1-4 |  cut -d "_" -f6- | cut -d "." -f2- )

# merging outG with sample.vcf
bcftools merge -r ${chr_pos} -m snps ${OUTG} ${chr} -Oz -o ${OUTDIR}/${SMP}_${chr}_anc.vcf.gz
bcftools index -t ${OUTDIR}/${SMP}_${chr}_anc.vcf.gz


