################################################################################
# Run RaxML on the T. forsythia SNP alignment for more robust tree building
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
        "tf.raxml.done",
        "tf.raxml.bs.done",
        "tf.raxml.support.done"

rule raxml:
    output:
        touch("tf.raxml.done")
    message: "Run raxml on tf SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/subset_Tf_snpAlignment.fasta"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --prefix tf \
         --threads {threads} \
         --seed 2
        """

rule bootstrap:
    input: 
        "tf.raxml.done"
    output:
        touch("tf.raxml.bs.done")
    message: "Run bootstrapping on trees for tf SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/subset_Tf_snpAlignment.fasta"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --bootstrap \
         --bs-trees 200 \
         --prefix tf \
         --threads {threads} \
         --seed 2
        """

# next step for branch support
rule branch_support:
    input: 
        "tf.raxml.bs.done"
    output:
        touch("tf.raxml.support.done")
    message: "Get branch support for RaxML tree"
    conda: "ENVS_raxml.yaml"
    params: 
        tree = "tf.raxml.bestTree",
        boots = "tf.raxml.bootstraps"
    threads: 2
    shell:
        """
        raxml-ng --support \
        --tree {params.tree} \
        --bs-trees {params.boots} \
        --prefix tf \
        --threads {threads} 
        """

