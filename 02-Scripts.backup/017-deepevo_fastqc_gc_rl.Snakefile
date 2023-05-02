################################################################################
# Run FASTQC to get read length and GC content
#
# Irina Velsko, 05/31/2022
################################################################################

from glob import glob
import os
import re

# workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/sample_gc_rl/deep_evo"
workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/sample_gc_rl/mann2018_world"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/mann2018_world/eager2_out/samtools/filter/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        expand("{sample}.human_filtered.zip", sample=SAMPLES.keys())

rule fastqc:
    output:
        "{sample}.human_filtered.zip"
    message: "Run FASTQC on {wildcards.sample}"
    group: "infoseq"
    params:
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        fastqc {params.fasta} -o /mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/sample_gc_rl/mann2018_world/FastQC
        """
