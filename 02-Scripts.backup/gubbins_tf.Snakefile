################################################################################
# Run gubbins on the Tannerella forsythia full alignment to look for regions of recombination
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
        "tf.filtered_polymorphic_sites.fasta",
        "tf.gubbins.masked.aln"

rule gubbins:
    output:
        "tf.filtered_polymorphic_sites.fasta"
    message: "Run gubbins on Tf full alignment"
    conda: "ENVS_gubbins.yaml"
    params: 
        infa = "tf_fullAlignment.fasta"
    threads: 8
    shell:
        """
         run_gubbins.py --prefix tf \
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
        "tf.gubbins.masked.aln"
    message: "Mask recombinant regions of Tf alignment"
    conda: "ENVS_gubbins.yaml"
    params: 
        infa = "tf_fullAlignment.fasta"
    shell:
        """
        mask_gubbins_aln.py --aln {params.infa} \
        --gff tf.recombination_predictions.gff \
        --out {output}
        """

