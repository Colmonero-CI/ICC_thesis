#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 0-12:00:00
#SBATCH -J VEP
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools bcftools vep

# chromosome name
chr=$1

# path to chromosome vcf file directory
SRCH=$2

# path to output directory
OUTDIR=$3

# find the chromosome vcf file you want to annotate
FILE=$(find ${SRCH} -maxdepth 1 -name "*${chr}*.gz"  -exec basename {} \;) 

# run vep on the chromosome vcf file
## here we are using the macaque genome as the reference annotated genome (Ensembl); change accordingly
## the cache directory is where the VEP reference annotated genome is stored; change according to your setup
## the --pick flag is used to select the most severe consequence for each variant; change according to your needs
## the --fork flag is used to thread VEP; change according to your setup
## the --vcf flag is used to output the annotated vcf file; change according to your needs; all my filtering scripts are based on vcf file output
    # but I adivise using the standard output
    
vep -i ${SRCH}/${FILE} --species macaca_mulatta --cache --dir $VEP_CACHE \
    --pick \
    -o ${OUTDIR}/${chr}_vep.ann.vcf.gz --vcf --compress_output bgzip \
    --fork 4

# index the output vcf
bcftools index -t ${OUTDIR}/${chr}_vep.ann.vcf.gz
