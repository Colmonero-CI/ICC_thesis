# Beta-psmc analysis

# Author: Colmonero-Costeira, I. based on Jansen, A.

This directory contains all the scripts that I've used to run beta-psmc.

- For this to work, you first need to install beta-psmc:
https://github.com/ChenHuaLab/Beta-PSMC (either in the home directory as I have or change the path accordingly in the scripts)

- These scripts assumes that the vcf files are divided by chromosomes, you need a file that lists the paths of these vcf files as specified in the vcf2fastq.sh script.  This is called chr_list.txt in my scripts.

- You need a file that lists all the chromosomes to analyse and their full range, three columns: chromosome\tstart\tend (i.e. chromosome,1,chromosome length). This is called chr_intervals.txt in my scripts.

- It's a sample-specific analysis, so I've written it to take a text file with the list of samples as command line argument. Then the sample IDs are parsed into the next steps of the pipeline. This is called sample_list.txt.

- The parse_psmc.py script is just a parser to put all the output in a plottable table, separating the different bootstraps runs etc.

- Remember to change the paths to all the beta-psmc executables.

Scripts in pipeline:

- launch_beta_psmc.sh: script to run vcf2fastq.sh across all samples in the sample_list.txt 
- vcf2fastq.sh: script to convert the vcf to a fasta file by sample, by chromosome; and launch beta_psmc.sh script
- beta_psmc.sh: script that converts the fasta files by chromosome into the correct psmc format (.psmcfa); and launches the bootstraping  (beta_bootstraps.sh), psmc and beta-psmc analyses for each sample; also parses the outputs (parse_pscmc.py by Jansen, A.)
- beta_bootstraps.sh runs bootstraping runs and parses the outputs (parse_pscmc.py by Jansen, A.)

To launch the pipeline:

- Keep the inputs in the launching directory; or change the scripts to the new directory architecture you want to use
- Launch as $sbatch launch_beta_psmc.sh PATH_TO_SAMPLE_LIST; e.g., $sbatch launch_beta_psmc.sh /proj/sllstore2017021/nobackup/IVO/betapsmc_ICC/sample_list.txt





