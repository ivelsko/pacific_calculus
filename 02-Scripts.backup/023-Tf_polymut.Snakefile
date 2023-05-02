################################################################################
# Polymorphic rate over protein-coding genes with polymut.py for Tannerella forsythia
#
# Irina Velsko 16/02/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/T_forsythia"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/T_forsythia/masked_mapped_bams/*.bam"):
	SAMPLES[os.path.basename(sample).split(".m")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        expand("output_masked_bams/{sample}.polymut.tsv", sample=SAMPLES.keys())
        
rule polymut:
    output:
        "output_masked_bams/{sample}.polymut.tsv"
    message: "Run polymut.py on {wildcards.sample} for Tf"
    params:
        refgff = "/mnt/archgen/microbiome_sciences/reference_genomes/Tannerella_forsythia/Tannerella_forsythia_9212_proka_format.gff",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate /mnt/archgen/users/huebner/automatic_MAG_refinement/conda/c7f0391e54d836540d3e1f91a1703ac7
        set -u
        
        polymut.py {params.bam} \
        --sortindex \
        --mincov 5 \
        --minqual 30 \
        --dominant_frq_thrsh 0.8 \
        --gff_file {params.refgff} > {output}
        """

