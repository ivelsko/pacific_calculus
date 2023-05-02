################################################################################
# Run SourceTracker on a rarefied species table for the Pacific data
# so it's consistent with the comparative data, which was all analyzed 
# at the species level
#
# Irina Velsko 22/07/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/sourcetracker/shotgun/species"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        "sourcetracker.done"
        
rule sourcetracker:
    output:
        "sourcetracker.done"
    message: "Run SourceTracker for the Pacific data at species level"
    params:
        otu = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/oceania_sourcetracker_table_species.txt",
        map = "/mnt/archgen/microbiome_calculus/pacific_calculus/05-Documentation.backup/source_tracker_mappingfile_20220722.tsv"
    shell:
        """
        Rscript \
        /projects1/users/velsko/bin/sourcetracker-1.0.1/sourcetracker_for_qiime.r \
        -i {params.otu} \
        -m {params.map} \
        -o shotgun_sourcetracker_species \
        -r 10000 \
        --train_rarefaction 5000 \
        -v
        """

