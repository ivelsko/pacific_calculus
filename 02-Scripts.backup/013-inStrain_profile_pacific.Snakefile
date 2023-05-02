################################################################################
# Prepare  Pacific calculus data to run inStrain
#
# Irina Velsko, 08/10/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/masked_mapped_bams/*.bam"):
	SAMPLES[os.path.basename(sample).split(".p")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("output/{sample}.IS", sample=SAMPLES.keys())

rule profile:
    output:
        directory("output/{sample}.IS")
    group: "dRep"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.fasta",
        bam = lambda wildcards: SAMPLES[wildcards.sample],
        genes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.genes.fna",
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.stb"
    message: "inStrain profile"
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate inStrain
        set -u

        inStrain profile {params.bam} {params.reffa} \
        -o {output} \
        -p {threads} \
        -g {params.genes} \
        -s {params.stb} \
        --min_insert 12 \
        --database_mode
        """
