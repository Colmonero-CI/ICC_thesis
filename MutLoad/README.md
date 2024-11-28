# Mutational load analysis

# Author: Colmonero-Costeira, I. 

This directory contains all the scripts that I've used to run mutational load analysis.

- This pipeline is devided in two sections. Polarization of the ancestral alleles and annotation of variants. 

- These scripts assumes that the vcf files are divided by chromosomes, you need : 
    - a file that lists the paths of these vcf files as specified in the out_to_ref.sh script. This is called chr_path.txt in my scripts.
    - a file that has the chromosome names as specified in the launch_BEP_ann.sh script. This is called chr_list.txt in my scripts.

## Polarization of the ancestral alleles

- For this to work, you first need to install minimap2:
https://github.com/lh3/minimap2 (either in the home directory as I have or change the path accordingly in the scripts).

- Other programs we will be using: samtools and bcftools.

- You need the .fna files for the reference genome (here it's GCF_003339765.1_Mmul_10) and for the outgroups (GCF_000951035.1_C_angolensis and GCA_028645565.1_P_papio). 
    -  You need a file with the paths for the outgroup genomes as specified in the out_to_ref.sh script. This is called outg_path.txt in my scripts.

- This pipeline assumes less then 20% divergence between the reference genome and the outgroups.

- Remember to change the paths to all the executables.

Scripts in the sub-pipeline:

- outg_to_ref: script to align the outgroups to the reference genome; and launch out_merge.sh for each chromosome vcf file
- outg_merge.sh: script to merge the outgroups to the population level vcf files.

To launch the sub-pipeline:

- Keep the inputs in the launching directory; or change the scripts to the new directory architecture you want to use
- Launch as $sbatch out_to_ref.sh PATH_TO_REF OUTPUT_PATH; e.g., $sbatch out_to_ref.sh /proj/sllstore2017021/nobackup/IVO/REF_GENOMES/GCF_003339765.1_Mmul_10/GCF_003339765.1_Mmul_10_genomic.fna /proj/sllstore2017021/nobackup/IVO/outg_align_merge_ICC

## Annotation of variants with VEP

- For this to work, you first need to install VEP: https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html.

- Other programs we will be using: bcftools.

- This sub-pipeline is designed for mutational load purposes (--pick flag in VEP_ann.sh script) but can be changed according to your annotation needs. See the VEP page for additional flags (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html).

- Remember to change the paths to all the executables.

Scripts in the sub-pipeline:

- launch_VEP_ann.sh: script to cript to run VEP_ann.sh across all chromosomes in chr_list.txt
- VEP_ann.sh: script to annotate the chromosome .vcf files using Ensembl's VEP

To launch the sub-pipeline:

- Keep the inputs in the launching directory; or change the scripts to the new directory architecture you want to use
- Launch as $sbatch launch_VEP_ann.sh PATH_TO_INPUT_VCFs OUTPUT_PATH; e.g., $sbatch launch_VEP_ann.sh /proj/sllstore2017021/nobackup/IVO/REF_GENOMES/Mmul_10_align/anc_call /proj/sllstore2017021/nobackup/IVO/REF_GENOMES/Mmul_10_align/anc_call/VEP_ann/GB/SEVERE





