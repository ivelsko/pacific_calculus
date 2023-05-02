################################################################################
# Test different substitution models for building the Abot439 tree from Pacific data
#
# Irina Velsko 28/03/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/modeltest_ng"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
         "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/modeltest_ng/modeltest_ng.out"
        
rule modeltest_ng:
    output:
        "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/modeltest_ng/modeltest_ng.out"
    message: "Infer the best-fitting substitution model using Modeltest-NG"
    conda: "ENVS_modeltest.yaml"
    resources:
        mem = 32,
        cores = 16
    params:
        fasta = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/subset_abot439_snpAlignment.fasta",
        prefix = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/modeltest_ng/modeltest_ng"
    threads: 16
    shell:
        """
        modeltest-ng -i {params.fasta} \
            --datatype nt \
            --output {params.prefix} \
            -p {threads} \
            -r 0 \
            -s 11
        """
