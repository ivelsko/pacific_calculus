################################################################################
# Run RaxML on the Abot439 SNP alignment for more robust tree building
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
        "abot439.raxml.done",
        "abot439.raxml.bs.done",
        "abot439.raxml.support.done"

rule raxml:
    output:
        touch("abot439.raxml.done")
    message: "Run raxml on Abot439 SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/subset_abot439_snpAlignment.fasta"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --prefix abot439 \
         --threads {threads} \
         --seed 2
        """

rule bootstrap:
    input:
        "abot439.raxml.done"
    output:
        touch("abot439.raxml.bs.done")
    message: "Run raxml on Abot439 SNP alignment"
    conda: "ENVS_raxml.yaml"
    params: 
        infa = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/subset_abot439_snpAlignment.fasta"
    threads: 8
    shell:
        """
         raxml-ng --msa {params.infa} \
         --model TN93+G4m \
         --bootstrap \
         --bs-trees 200 \
         --prefix abot439 \
         --threads {threads} \
         --seed 2
        """


# next step for branch support
rule branch_support:
    input: 
        "abot439.raxml.bs.done"
    output:
        touch("abot439.raxml.support.done")
    message: "Get branch support for RaxML tree"
    conda: "ENVS_raxml.yaml"
    params: 
        tree = "abot439.raxml.bestTree",
        boots = "abot439.raxml.bootstraps"
    threads: 2
    shell:
        """
        raxml-ng --support \
        --tree {params.tree} \
        --bs-trees {params.boots} \
        --prefix abot439 \
        --threads {threads} 
        """

