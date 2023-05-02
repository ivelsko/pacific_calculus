#!/usr/bin/env bash

sbatch \
-c 4 \
--mem 32000 \
-o ~/slurm_logs/slurm.%j.out \
-e ~/slurm_logs/slurm.%j.err \
-p medium \
--mail-type=fail \
--mail-user=fagernaes@shh.mpg.de \
--wrap="/projects1/tools/java/jdk-11.0.2/bin/java \
-jar /projects1/clusterhomes/huebler/RMASifter/AMPS/MaltExtract1.7.jar \
-input /projects1/microbiome_calculus/pacific_calculus/04-Analysis/malt/output_fullNT2017/ \
-output /projects1/microbiome_calculus/pacific_calculus/04-Analysis/eukaryotes/ \
-m me_po \
-c maltex_config.txt"
