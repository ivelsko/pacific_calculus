################################################################################
# Prepare  Pacific calculus data to run inStrain
#
# Irina Velsko, 08/10/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/masked_mapped_bams/*.bam"):
	SAMPLES[os.path.basename(sample).split(".p")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("inStrain_output/{sample}.is12.IS", sample=SAMPLES.keys()),
        expand("inStrain_output/{sample}.is24.IS", sample=SAMPLES.keys()),
        expand("inStrain_output/{sample}.is36.IS", sample=SAMPLES.keys()),
        expand("inStrain_output/{sample}.is48.IS", sample=SAMPLES.keys())

rule profile_is12:
    output:
        directory("inStrain_output/{sample}.is12.IS")
    group: "dRep"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
        bam = lambda wildcards: SAMPLES[wildcards.sample],
        genes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fna",
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.stb"
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

rule profile_is24:
    output:
        directory("inStrain_output/{sample}.is24.IS")
    group: "dRep"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
        bam = lambda wildcards: SAMPLES[wildcards.sample],
        genes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fna",
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.stb"
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
        --min_insert 24 \
        --database_mode
        """

rule profile_is36:
    output:
        directory("inStrain_output/{sample}.is36.IS")
    group: "dRep"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
        bam = lambda wildcards: SAMPLES[wildcards.sample],
        genes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fna",
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.stb"
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
        --min_insert 36 \
        --database_mode
        """

rule profile_is48:
    output:
        directory("inStrain_output/{sample}.is48.IS")
    group: "dRep"
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
        bam = lambda wildcards: SAMPLES[wildcards.sample],
        genes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fna",
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.stb"
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
        --min_insert 48 \
        --database_mode
        """
