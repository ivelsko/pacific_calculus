################################################################################
# Map Pacific calculus against Abot439 reference genome for polymut,
# and perform coverage calculations
#
# Irina Velsko 04/11/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut" 

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/screening_all_rerun/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("bwa_mapping/{sample}.abot439.m_s_r.cov", sample=SAMPLES.keys()),
        expand("bwa_mapping/{sample}.abot439.m_s_r.bam", sample=SAMPLES.keys()),
        expand("bwa_mapping/{sample}.abot439.m_s_r.bam.bai", sample=SAMPLES.keys()),
        expand("masked_mapped_bams/{sample}.abot439.m_s_r.mask_1_1.bam", sample=SAMPLES.keys()),
        expand("masked_mapped_bams/{sample}.abot439.m_s_r.mask_13_13.bam", sample=SAMPLES.keys())

rule bwa_aln:
    output:
        temp("bwa_mapping/{sample}.abot439.sai")
    message: "Align {wildcards.sample} against Abot439 reference genome using BWA aln"
    params: 
        reffa = "/mnt/archgen/microbiome_sciences/reference_genomes/Anaerolineaceae_b_oral_taxon_439/Anaerolineaceae_b_oral_taxon_439.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 4
    shell:
        """
        bwa aln -n 0.01 -l 1024 -t {threads} \
            {params.reffa} \
            {params.fastq} > {output}
        """

rule bwa_samse:
    input:
        "bwa_mapping/{sample}.abot439.sai"
    output:
        "bwa_mapping/{sample}.abot439.m_s.bam"
    message: "Generate alignment file for {wildcards.sample} against Abot439 referene genome"
    params:
        reffa = "/mnt/archgen/microbiome_sciences/reference_genomes/Anaerolineaceae_b_oral_taxon_439/Anaerolineaceae_b_oral_taxon_439.fna",
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bwa samse \
            {params.reffa} \
            {input} \
            {params.fastq} | \
        samtools view -Sb -F 4 - | samtools sort - -o {output} 
        """

rule samtools_rmdup:
    input:
        "bwa_mapping/{sample}.abot439.m_s.bam"
    output:
        "bwa_mapping/{sample}.abot439.m_s_r.bam"
    message: "Remove duplicate mapped reads for {wildcards.sample} against Abot439 referene genome"
    params:
    shell:
        """
        samtools rmdup -s {input} {output}
        """

rule samtools_index:
    input:
        "bwa_mapping/{sample}.abot439.m_s_r.bam"
    output:
        "bwa_mapping/{sample}.abot439.m_s_r.bam.bai"
    message: "Index bam file for {wildcards.sample} mapped against Abot439 referene genome"
    params:
    shell:
        """
        samtools index {input}
        """

rule bedtools_coverage:
    input:
        "bwa_mapping/{sample}.abot439.m_s_r.bam"
    output:
        "bwa_mapping/{sample}.abot439.m_s_r.cov"
    message: "Calculate depth/breadth of coverage of {wildcards.sample} mapped against Abot439 referene genome"
    params:
        refbed = "/mnt/archgen/microbiome_calculus/Abot439/04-analysis/reference_genome_mapping/abot439_reference.bed"
    shell:
        """
        bedtools coverage -a {params.refbed} -b {input} -hist > {output}
        """

rule mask_1:
    input:
        "bwa_mapping/{sample}.abot439.m_s_r.bam"
    output:
        "masked_mapped_bams/{sample}.abot439.m_s_r.mask_1_1.bam"
    message: "Mask ends of mapped udg-half {wildcards.sample}"
    params:
    shell:
        """
        python /mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/maskTerminalDeam.py \
        -i {input} \
        -l 1 \
        -r 1 \
        -o {output}
        """

rule mask_13:
    input:
        "bwa_mapping/{sample}.abot439.m_s_r.bam"
    output:
        "masked_mapped_bams/{sample}.abot439.m_s_r.mask_13_13.bam"
    message: "Mask ends of mapped udg-half {wildcards.sample}"
    params:
    shell:
        """
        python /mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/maskTerminalDeam.py \
        -i {input} \
        -l 13 \
        -r 13 \
        -o {output}
        """
