################################################################################
# Run gubbins on the Abot439 full alignment to look for regions of recombination
# (to remove from the SNP aligment-based trees?)
#
# Irina Velsko 18/06/2024
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/gubbins"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

rule all:
    input: 
        "abot439.filtered_polymorphic_sites.fasta",
        "abot439.gubbins.masked.aln"

rule gubbins:
    output:
        "abot439.filtered_polymorphic_sites.fasta"
    message: "Run gubbins on Abot439 full alignment"
    conda: "ENVS_gubbins.yaml"
    params: 
        infa = "abot439_fullAlignment.fasta"
    threads: 8
    shell:
        """
         run_gubbins.py --prefix abot439 \
         --first-tree-builder fasttree \
         --tree-builder raxml \
         --best-model \
         --bootstrap 100 \
         --threads {threads} \
         --verbose \
         {params.infa}
        """


rule mask_alignment:
    output:
        "abot439.gubbins.masked.aln"
    message: "Mask recombinant regions of Abot439 alignment"
    conda: "ENVS_gubbins.yaml"
    params: 
        infa = "abot439_fullAlignment.fasta"
    shell:
        """
        mask_gubbins_aln.py --aln {params.infa} \
        --gff abot439.recombination_predictions.gff \
        --out {output}
        """

