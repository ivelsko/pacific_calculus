################################################################################
# Map gargammel-simulated UDG-half reads against each ZymoBIOMICS reference genome 
#
# Irina Velsko, 05/31/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/eager_output/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.bam.bai", sample=SAMPLES.keys()),
        expand("zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.bam", sample=SAMPLES.keys())

rule aln:
    output:
        temp("zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.sai")
    params:
        refgenome = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome.fasta",
        infile = lambda wildcards: SAMPLES[wildcards.sample]
    message: "bwa map {wildcards.sample} against ZymoBIOMICS Cryptococcus neoformans reference genome"
    shell:
        """
        bwa aln \
        -n 0.02 \
        -l 1024 \
        {params.refgenome} \
        {params.infile} > {output}
        """

rule samse:
    input: 
        "zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.sai"
    output:
        "zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.bam"
    params:
        refgenome = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/ZymoBIOMICS.STD.refseq.v2/Genomes/Cryptococcus_neoformans_draft_genome.fasta",
        infile = lambda wildcards: SAMPLES[wildcards.sample]
    message: "bwa samse {wildcards.sample} against ZymoBIOMICS Cryptococcus neoformans reference genome"
    shell:
        """
        bwa samse \
        {params.refgenome} \
        {input} \
        {params.infile}  | samtools view -Sb -F 4 - | samtools sort - -o {output}
        """

rule index:
    input:
        "zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.bam"
    output:
        "zymo_individual_mapping/{sample}.Cryptococcus_neoformans_draft_genome.bam.bai"
    params:
    message: "index {wildcards.sample} against ZymoBIOMICS Cryptococcus neoformans reference genome"
    shell:
        """
        samtools index {input}
        """
