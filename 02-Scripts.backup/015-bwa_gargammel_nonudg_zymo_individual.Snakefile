################################################################################
# Map gargammel-simulated non-UDG reads against each ZymoBIOMICS reference genome 
#
# Irina Velsko, 05/31/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/ZymoBIOMICS.STD.refseq.v2/Genomes/*.fasta"):
	SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("zymo_individual_mapping/non_udg.{sample}.bam", sample=SAMPLES.keys()),
        expand("zymo_individual_mapping/non_udg.{sample}.bam.bai", sample=SAMPLES.keys())

rule aln:
    output:
        temp("zymo_individual_mapping/non_udg.{sample}.sai")
    params:
        refgenome = lambda wildcards: SAMPLES[wildcards.sample],
        infile = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/eager_output/samtools/filter/non_udg.unmapped.fastq.gz"
    message: "bwa map non-udg gargammel against ZymoBIOMICS reference genome {wildcards.sample}"
    shell:
        """
        bwa aln \
        -n 0.02 \
        -l 1024 \
        {params.refgenome} \
        {params.infile} > {output}
        """

rule samse:
    input: "zymo_individual_mapping/non_udg.{sample}.sai"
    output:
        "zymo_individual_mapping/non_udg.{sample}.bam"
    params:
        refgenome = lambda wildcards: SAMPLES[wildcards.sample],
        infile = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/eager_output/samtools/filter/non_udg.unmapped.fastq.gz"
    message: "bwa samse non-udg gargammel against ZymoBIOMICS reference genome {wildcards.sample}"
    shell:
        """
        bwa samse \
        {params.refgenome} \
        {input} \
        {params.infile}  | samtools view -Sb -F 4 - | samtools sort - -o {output}
        """

rule index:
    input: 
        "zymo_individual_mapping/non_udg.{sample}.bam"
    output:
        "zymo_individual_mapping/non_udg.{sample}.bam.bai"
    params:
    message: "index non_udg mapped against {wildcards.sample}"
    shell:
        """
        samtools index {input}
        """
