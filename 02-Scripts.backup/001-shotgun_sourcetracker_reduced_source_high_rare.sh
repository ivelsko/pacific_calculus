#!/usr/bin/env bash

OUTDIR=/projects1/microbiome_calculus/pacific_calculus/04-Analysis/sourcetracker/shotgun

sbatch \
	-c 10 \
	--mem=24000 \
	-o ~/slurm_logs/slurm.%j.out \
	-e ~/slurm_logs/slurm.%j.err \
	--partition=long \
	--mail-type=fail \
	--mail-type=time_limit \
	--mail-user=fagernaes@shh.mpg.de \
	--export=ALL \
	-J "Sourcetracker" \
	--wrap="Rscript \
	/projects1/users/velsko/bin/sourcetracker-1.0.1/sourcetracker_for_qiime.r \
	-i /projects1/microbiome_calculus/pacific_calculus/05-Documentation.backup/otu_table_sourcetracker_genus_20210617.txt \
	-m /projects1/microbiome_calculus/pacific_calculus/05-Documentation.backup/oceania_mappingfile_20210611.tsv \
	-o "$OUTDIR"/shotgun_sourcetracker_genus_table_log_"$(date +"%Y%m%d")" \
	-r 10000 \
	--train_rarefaction 5000 \
	-v"

