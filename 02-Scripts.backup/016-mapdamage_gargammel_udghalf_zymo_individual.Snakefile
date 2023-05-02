################################################################################
# Run mapdamage on the gargammel udg-half sample mapped to each ZymoBIOMICS reference genome
#
# Irina Velsko, 05/31/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/mapdamage_gargammel_zymo_individual/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/ZymoBIOMICS.STD.refseq.v2/Genomes/*.fasta"):
	SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("udg_half_{sample}", sample=SAMPLES.keys())

rule mapdamage:
    output:
        directory("udg_half_{sample}")
    params:
        bam = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/zymo_individual_mapping/udg_half.{sample}.bam",
        reffa = lambda wildcards: SAMPLES[wildcards.sample]
    message: "mapdamage gargammel udg-half reads"
    threads: 12
    shell:
        """
        mapDamage -i {params.bam} \
        -r {params.reffa} \
        -d {output}
        """
