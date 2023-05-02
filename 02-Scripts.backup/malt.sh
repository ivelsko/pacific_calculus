#!/usr/bin/env bash

## Set your project parent directory and MALT-related directories
PROJDIR=/projects1/microbiome_calculus/pacific_calculus
MALT_RUN=/projects1/malt/versions/malt040/malt-run

sbatch \
-c 112 \
--mem=1850000 \
--partition=supercruncher \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
--mail-type=fail \
--mail-type=time_limit \
--mail-user=fagernaes@shh.mpg.de \
-J "pacific_all_MALT" \
--wrap="$PROJDIR/02-Scripts.backup/007-malt-genbank-nt_2017_2step_85pc_supp_0.01 $PROJDIR/03-Preprocessing/screening_all_rerun/samtools/filter/*unmapped.fastq.gz $PROJDIR/04-Analysis/malt/output_fullNT2017/"

