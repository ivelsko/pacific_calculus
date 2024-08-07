################################################################################
# Run dRep on Anaerolineaceae genomes generated by mapping Pacific 
# calculus against the reference genome
#
# Irina Velsko, 17/07/2024
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/multivcfanalyzer_hom/dRep"


if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        "out/data_tables/Ndb.csv"


rule dRep:
    output:
        "out/data_tables/Ndb.csv"
    message: "Run dRep to cluster mapping-based Anaerolineaceae genomes from the Pacific calculus project"
    params: 
        genomes = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/phylogenies/anaerolineaceae_bacterium_ot_439/multivcfanalyzer_hom/dRep/abot439_pacific_mapping_genome_list.tsv"
    threads: 8
    shell:
        """
        dRep dereplicate -p {threads} out/ -g {params.genomes} -comp 40 --S_algorithm ANImf -pa 0.95 -sa 0.99
        """
