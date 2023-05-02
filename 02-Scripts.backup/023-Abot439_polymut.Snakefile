################################################################################
# Polymorphic rate over protein-coding genes with polymut.py for Abot439
#
# Irina Velsko 01/12/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/masked_mapped_bams/*.bam"):
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
    message: "Run polymut.py on {wildcards.sample} for Abot439"
    params:
        refgff = "/mnt/archgen/microbiome_sciences/reference_genomes/Anaerolineaceae_b_oral_taxon_439/Anaerolineaceae_b_oral_taxon_439.gff",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate /home/alexander_huebner/miniconda3/envs/cmseq_Genbank
        set -u
        
        polymut.py {params.bam} \
        --sortindex \
        --mincov 5 \
        --minqual 30 \
        --dominant_frq_thrsh 0.8 \
        --gff_file {params.refgff} > {output}
        """

