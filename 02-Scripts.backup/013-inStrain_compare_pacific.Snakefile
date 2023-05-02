################################################################################
# Run inStrain compare on Pacific calculus data 
#
# Irina Velsko, 05/29/2022
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/output/*"):
	SAMPLES[os.path.basename(sample).split("z")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        directory("output/pac_reference_genomes.IS.COMPARE")

rule compare:
    input:
        expand("output/{sample}", sample=SAMPLES)
    output:
        directory("output/pac_reference_genomes.IS.COMPARE")
    group: "dRep"
    params:
        stb = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/inStrain/reference_fasta/pac_reference_genomes.stb"
    message: "inStrain compare"
    threads: 12
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate inStrain
        set -u

        inStrain compare -i {input} \
        -o {output} \
        -p {threads} \
        -s {params.stb} \
        --database_mode
        """

