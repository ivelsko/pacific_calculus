################################################################################
# Prepare gargammel simulated data for Pacific project to run inStrain by
# mapping to the cat'd ZymoBIOMICS reference genome file
#
# Irina Velsko, 05/31/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/zymo_all_mapping"

#### SAMPLES ###################################################################
SAMPLES = {os.path.basename(fn).replace("_1.fastq.gz", ""): fn.replace("_1.fastq.gz", "")
		   for fn in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/eager_output/fastqs/per_sample/*_1.fastq.gz")}
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.zymo_reference_genomes.bam", sample=SAMPLES.keys())


rule bwa_aln:
    output:
        temp("{sample}_{i}.sai")
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
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
        pipe("{sample}.zymo_reference_genomes.sam")
    params:
        reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/zymo_reference_genome/zymo_ref_genomes.fasta",
        reads = lambda wildcards: " ".join([f"{SAMPLES[wildcards.sample]}_{i}.fastq.gz" for i in range(1, 3)])
    message: "bwa sampe"
    shell:
        """
        bwa sampe {params.reffa} {input} {params.reads} > {output}
        """

rule bam:
    input:
        "{sample}.zymo_reference_genomes.sam"
    output:
        "{sample}.zymo_reference_genomes.bam"
    message: "sam2bam"
    shell:
        """
        samtools view -S -b {input} > {output}
        """
