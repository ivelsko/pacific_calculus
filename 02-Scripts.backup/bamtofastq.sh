#!/usr/bin/env bash

for FILE in /projects1/microbiome_calculus/pacific_calculus/01-Data/non_human/batch2/*/*extractunmapped.bam;
do sbatch -c 4 --mem 32000 --partition=short -J "convert" --wrap="samtools fastq $FILE >> $FILE.fastq && gzip $FILE.fastq"
done
