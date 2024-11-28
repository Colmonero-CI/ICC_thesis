#!/bin/bash

#SBATCH -A naiss2023-5-506
#SBATCH -p core
#SBATCH -n 1
#SBATCH --array=1-20
#SBATCH -t 0-10:00:00
#SBATCH -J vcf2fasta
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=costeira.ivo@gmail.com
#SBATCH -o ./logs/%x-%j.out
#SBATCH -e ./logs/%x-%j.error

module load bioinfo-tools bcftools samtools

# sample as command line input
SAMPLE=$1

# file listing the vcf files, one per line
vcf_files=chr_list.txt #change to the correct path

# get the vcf file for this array id
VCF=$(cat ${vcf_files} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# tempdir for the sample
mkdir -p ${SAMPLE}_infiles
# path to reference genome
REF=/proj/sllstore2017021/nobackup/IVO/REF_GENOMES/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna #change to the correct path

# to start we need to convert the vcf to a fasta file by sample, by chromosome
# get the chromosome from the vcf
CHROM=$(bcftools view -H ${VCF} | head -n 1 | cut -f 1 )
#set start 
START=$(cat chrm_auto_intervals.txt | awk -v chrom=${CHROM} '$1 == chrom {print $2}')

# get length of the chromosome
END=$(cat chrm_auto_intervals.txt | awk -v chrom=${CHROM} '$1 == chrom {print $3}')

echo ${CHROM}:${START}-${END} > ${SAMPLE}_infiles/FA_${CHROM}.txt

FA=/proj/sllstore2017021/nobackup/IVO/betapsmc_ICC/${SAMPLE}_infiles/FA_${CHROM}.txt

# subset the reference genome by chromosome and convert this samples genotype to a fasta file
samtools faidx ${REF} -r ${FA} | \
bcftools consensus -s ${SAMPLE} --haplotype I ${VCF} --absent "N" --missing "N"  | gzip > ${SAMPLE}_infiles/${CHROM}_${SAMPLE}.fa.gz

# when all the jobs are done, submit the psmc job 
if [ ${SLURM_ARRAY_TASK_ID} -eq 20 ]
then
	sbatch --dependency=afterok:${SLURM_JOB_ID} beta_psmc.sh ${SAMPLE}
fi