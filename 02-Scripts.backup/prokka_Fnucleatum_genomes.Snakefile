################################################################################
# Run prokka on all F. nucleatum genomes to use as input for polymut
#
# Irina Velsko 23/03/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/*.fna"):
	SAMPLES[os.path.basename(sample).split(".f")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")

rule all:
    input:
        expand("prokka/{sample}.gff.gz", sample=SAMPLES.keys())


rule prokka:
    output:
        "prokka/{sample}.gff.gz"
    message: "Run Prokka on contigs: {wildcards.sample}"
    resources:
        mem = 16,
        cores = 8
    params:
        outdir = "prokka",
        fasta = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        prokka --outdir {params.outdir} \
               --prefix {wildcards.sample} \
               --force \
               --compliant \
               --metagenome \
               --cpus {threads} \
               --debug \
               {params.fasta}
        pigz -f -p 4 {params.outdir}/{wildcards.sample}.{{faa,ffn,fna,gbk,gff,tsv,txt}}
        """
