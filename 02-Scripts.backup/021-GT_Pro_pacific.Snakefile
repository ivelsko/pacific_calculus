################################################################################
# Prepare  Pacific calculus data to run inStrain
#
# Irina Velsko, 08/10/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/GT_Pro/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/screening_all_rerun/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".p")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.IS", sample=SAMPLES.keys())

rule gt_pro_genotype:
    output:
        directory("{sample}.IS")
    group: "dRep"
    params:
        db = "/mnt/archgen/microbiome_sciences/reference_databases/built/GT_Pro",
        infile = lambda wildcards: SAMPLES[wildcards.sample],
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.stb"
    message: "inStrain profile"
    threads: 12
    shell:
        """
        /mnt/archgen/users/velsko/bin/gt-pro/GT_Pro genotype \
        -d {db} \
        -t {threads} \
        {infile}
        """
