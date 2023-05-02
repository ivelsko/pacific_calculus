################################################################################
# Mask the ends of mapped reads to "hide" aDNA damage patterns from SNP 
# detecting programs on Pacific calculus data, using Alex's script
#
# Irina Velsko, 17/06/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/masked_mapped_bams"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/bwa_mapping/*.bam"):
	SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("{sample}.mask_1_1.bam", sample=SAMPLES.keys()),
        expand("{sample}.mask_13_13.bam", sample=SAMPLES.keys())

rule mask1:
    output:
        "{sample}.mask_1_1.bam"
    params:
        inbam = lambda wildcards: SAMPLES[wildcards.sample]
    message: "mask ends of mapped reads for udg-half {wildcards.sample}"
    shell:
        """
        python /mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/maskTerminalDeam.py \
        -i {params.inbam} \
        -l 1 \
        -r 1 \
        -o {output}
        """

rule mask13:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        "{sample}.mask_13_13.bam"
    params:
    message: "mask ends of mapped reads for non-udg {wildcards.sample}"
    shell:
        """
        python /mnt/genotyping/sk_pipelines/projects/aDNA_Flores/scripts/maskTerminalDeam.py \
        -i {input} \
        -l 13 \
        -r 13 \
        -o {output}
        """
