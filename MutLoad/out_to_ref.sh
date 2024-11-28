#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p node
#SBATCH -N 1
#SBATCH -t 5:00:00
#SBATCH -J minimap2
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools minimap2 samtools bcftools

# path to reference
REF=$1

# path to output directory
OUTDIR=$2

# output name
SMP=$3

# txt file with outgroup paths
OUTG=outg_path.txt

# txt with chromosome vcf files paths
CHR=chr_path.txt

for out in ${OUTG}
do
    # align outgroup to reference using minimap2
    ## option asm20 assumes less that 20% divergence between the reference genome and the outgroup.
    ## output is a sam file that is sorted and converted to bam format
    minimap2 -ax asm20 ${REF} ${out} | samtools sort | samtools view -bh -o ${SNIC_TMP}/${out}_mmul10.bam
    samtools index ${SNIC_TMP}/${out}_mmul10.bam

    # call all variants on the outgroups using bcftools
    bcftools mpileup -f ${REF} ${SNIC_TMP}/${out}_mmul10.bam | bcftools call -m -Oz -o ${OUTDIR}/${out}_mmul10.vcf.gz
    bcftools index -t ${OUTDIR}/${out}_mmul10.vcf.gz

    # same as above but keep only keep genotypes candidate for ancestral allele states (0/0 genotypes) 
    ## use only if doing mutational load analyses
    #bcftools mpileup -f ${REF} ${SNIC_TMP}/${out}_mmul10.bam | bcftools call -m | bcftools view -i 'GT="0/0"' -Oz -o ${OUTDIR}/${out}_anc.vcf.gz
    #bcftools index -t ${OUTDIR}/${out}_anc.vcf.gz
done

# merge outgroups
# create a list of outgroup vcf files
find . -maxdepth 1 -name "*.vcf.gz"  -exec basename {} \; | sort -V > outg_to_ref_list.txt 

# merge outgroups, keep only snps and remove missing genotypes
bcftools merge -m snps -l outg_to_ref_list.txt | bcftools view -e 'GT[*]="mis"' -Oz -o ${OUTDIR}/${SMP}.vcf.gz
bcftools index -t ${OUTDIR}/${SMP}.vcf.gz

# launch merging of outgroup and population level vcf files (by chromosome)

for chr in ${CHR}
do 
    sbatch out_merge.sh ${chr} ${OUTG} ${SRCH} ${OUTDIR} ${SMP}
done