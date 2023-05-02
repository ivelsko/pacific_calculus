################################################################################
# Prepare simulated data using gargammel for testing inStrain on ancient reads
# UDG-half-treated reads
#
# Irina Velsko, 05/30/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel"

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        directory("output/udg_half")

rule gargammel:
    output:
        directory("output/udg_half")
    params:
        infolder = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/input/",
        misinc = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/mapdamage_profiles_for_gargammel/results_EFE002.B0101.mapped/dnacomp.txt",
        mapdamage = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/mapdamage_profiles_for_gargammel/results_EFE002.B0101.mapped/misincorporation.txt",
        fragments = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/EFE002.B0101_lgdist.txt",
        ss = "HS25"
    message: "gargammel udg-half"
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate gargammel
        set -u

        gargammel \
        --comp 0.99,0,0.01 \
        -n 10000000 \
        --misince {params.misinc} \
        --misincb {params.misinc} \
        -f {params.fragments} \
        -mapdamagee {params.mapdamage} single \
        -mapdamageb {params.mapdamage} single \
        -rl 75 \
        -ss {params.ss} \
        -o {output} \
        {params.infolder}
        """

#         unset PERL5LIB
#         /mnt/archgen/users/velsko/bin/gargammel/gargammel.pl \
