################################################################################
# Get ANI of F. nucleatum subspecies genomes  with dRep
#
# Irina Velsko 03/04/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/F_nucleatum_gargammel/dRep"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

# leave the input with the wrong path so it doesn't delete the file when it "fails"
rule all:
    input:
        "dRep_out/data_table/Ndb.csv"

rule compare:
    output:
        "dRep_out/data_table/Ndb.csv"
    message: "cluster Fn subspecies genomes with dRep"
    params:
        bins = "Fn_genome_list.tsv",
    threads: 12
    shell:
        """
        dRep compare dRep_out/ -g {params.bins} --S_algorithm ANImf -pa 0.95 -sa 0.99
        """
