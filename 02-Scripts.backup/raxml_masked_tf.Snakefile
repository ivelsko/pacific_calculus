################################################################################
# Run RaxML on the masked Tannerella forsythia SNP alignment from gubbins
#
# Irina Velsko 18/06/2024
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/raxml"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

rule all:
    input: 
        "tf_masked.raxml.done",
        "tf_masked.raxml.bs.done",
        "tf_masked.raxml.support.done"

rule raxml:
    output:
        touch("tf_masked.raxml.done")
    message: "Run raxml on tf SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/gubbins/tf.gubbins.masked.aln"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --prefix tf_masked \
         --threads {threads} \
         --seed 2
        """

rule bootstrap:
    input: 
        "tf_masked.raxml.done"
    output:
        touch("tf_masked.raxml.bs.done")
    message: "Run raxml on tf SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/gubbins/tf.gubbins.masked.aln"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --bootstrap \
         --bs-trees 200 \
         --prefix tf_masked \
         --threads {threads} \
         --seed 2
        """

# next step for branch support
rule branch_support:
    input: 
        "tf_masked.raxml.bs.done"
    output:
        touch("tf_masked.raxml.support.done")
    message: "Get branch support for RaxML tree"
    conda: "ENVS_raxml.yaml"
    params: 
        tree = "tf_masked.raxml.bestTree",
        boots = "tf_masked.raxml.bootstraps"
    threads: 2
    shell:
        """
        raxml-ng --support \
        --tree {params.tree} \
        --bs-trees {params.boots} \
        --prefix tf_masked \
        --threads {threads} 
        """

