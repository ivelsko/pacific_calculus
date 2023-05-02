#!/usr/bin/env bash
 
#SBATCH -c 4                      # number of CPUs (here 4)
#SBATCH --mem=32G               # memory pool for all cores (here 32GB)
#SBATCH -o slurm.%j.out           # STDOUT (the standard output stream) into file slurm.<JOB_NUMBER>.out
#SBATCH -e slurm.%j.err           # STDERR (the output stream for errors) into file slurm.<JOB_NUMBER>.err
#SBATCH -p short                  # The queue or 'partition' you want to submit to, here the 'short' queue

java -Xmx16G -jar /projects1/tools/multivcfanalyzer/0.0.87/MultiVCFanalyzer_0-87.jar \
NA \
/projects1/microbiome_sciences/reference_genomes/Pseudopropionibacterium_propionicum/Pseudopropionibacterium_propionicum_F0230a.fa \
NA \
/projects1/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/pseudopropionibacterium_propionicum/multivcfanalyzer_hom/ \
T \
30 \
5 \
0.9 \
0.9 \
NA \
/projects1/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/pseudopropionibacterium_propionicum/genotyping/*/*.vcf.gz
