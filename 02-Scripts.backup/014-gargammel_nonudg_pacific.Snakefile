################################################################################
# Prepare simulated data using gargammel for testing inStrain on ancient reads
# Reads with no UDG treatment
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
        directory("output/non_udg")

rule gargammel:
    output:
        directory("output/non_udg")
    params:
        infolder = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/input/",
        misinc = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/mapdamage_profiles_for_gargammel/results_HCLVMBCX2-3505-07-00-01_S7.mapped/dnacomp.txt",
        mapdamage = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/mapdamage_profiles_for_gargammel/results_HCLVMBCX2-3505-07-00-01_S7.mapped/misincorporation.txt",
        fragments = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/benchmarking_ancient/gargammel/HCLVMBCX2-3505-07-00-01_S7_lgdist.txt",
        ss = "HS25"
    message: "gargammel non-udg"
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
        {params.infolder}/ 
        """
