################################################################################
# Prepare  Pacific calculus data to run inStrain
#
#
# Irina Velsko, 08/10/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/bwa_mapping"

#### SAMPLES ###################################################################
SAMPLES = {os.path.basename(fn).replace("_1.fastq.gz", ""): fn.replace("_1.fastq.gz", "")
		   for fn in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/input_fastq/*_1.fastq.gz")}
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.pac_reference_genomes.bam", sample=SAMPLES.keys()),
        "../reference_fasta/pac_reference_genomes.genes",
        expand("{sample}.pac_reference_genomes.bam.bai", sample=SAMPLES.keys()),
        expand("{sample}.pac_reference_genomes.cov", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("{sample}_{i}.sai")
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.fasta",
        reads = lambda wildcards: f"{SAMPLES[wildcards.sample]}_{wildcards.i}.fastq.gz"
    message: "bwa map R{wildcards.i}"
    shell:
        """
        bwa aln -n 0.01 -l 1024 {params.reffa} {params.reads} > {output}
        """

rule sampe:
    input:
        lambda wildcards: [f"{wildcards.sample}_{i}.sai" for i in range(1, 3)]
    output:
        pipe("{sample}.pac_reference_genomes.sam")
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.fasta",
        reads = lambda wildcards: " ".join([f"{SAMPLES[wildcards.sample]}_{i}.fastq.gz" for i in range(1, 3)])
    message: "bwa sampe"
    shell:
        """
        bwa sampe {params.reffa} {input} {params.reads} > {output}
        """

rule bam:
    input:
        "{sample}.pac_reference_genomes.sam"
    output:
        "{sample}.pac_reference_genomes.bam"
    message: "sam2bam"
    shell:
        """
        samtools view -S -b {input} | samtools sort - -o {output} 
        """

rule samtools_index:
    input:
        "{sample}.pac_reference_genomes.bam"
    output:
        "{sample}.pac_reference_genomes.bam.bai"
    message: "Index bam file for {wildcards.sample}"
    params:
    shell:
        """
        samtools index {input}
        """

rule prodigal:
    input:
        "../reference_fasta/pac_reference_genomes.fasta"
    output:
        "../reference_fasta/pac_reference_genomes.genes"
    params:
        faa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.proteins.faa",
        fna = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.genes.fna"
    message: "run prodigal"
    shell:
        """
        prodigal -i {input} -o {output} -a {params.faa} -d {params.fna} -p meta
        """

rule bedtools_coverage:
    input:
        "{sample}.pac_reference_genomes.bam"
    output:
        "{sample}.pac_reference_genomes.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against each genome"
    params:
        refbed = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.bed",
    shell:
        """
        bedtools coverage -a {params.refbed} -b {input} -hist > {output}
        """
