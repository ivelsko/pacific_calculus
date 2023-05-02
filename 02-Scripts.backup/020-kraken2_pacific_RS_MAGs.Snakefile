################################################################################
# Run CMC and comparative data through Kraken2 with 2 databases - 
# one RefSeq genomes only, one RefSeq genomes + MAGs from Pasolli, et al. 2019
#
# Irina Velsko 23/03/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/kraken2"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/kraken2/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("output/{sample}.kraken2.output.rspm.tsv", sample=SAMPLES.keys()),
        expand("output/{sample}.kraken2.report_mpa.rspm.tsv", sample=SAMPLES.keys()),
        expand("output/{sample}.kraken2.output.rs.tsv", sample=SAMPLES.keys()),
        expand("output/{sample}.kraken2.report_mpa.rs.tsv", sample=SAMPLES.keys())


rule kraken_rspm:
    output:
        outfmt = "output/{sample}.kraken2.output.rspm.tsv",
        repfmt = "output/{sample}.kraken2.report_mpa.rspm.tsv"
    message: "Run {wildcards.sample} through kraken2 database w/bacteria, archaea, human"
    params: 
        database = "/mnt/archgen/microbiome_sciences/reference_databases/built/refseq20191017_Pasolli2019/kraken2_db/MiniKraken_RefSeq1910PlusPasolliSGBs/",
        infile = lambda wildcards: SAMPLES[wildcards.sample],
    threads: 32
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate kraken2
        set -u

        kraken2 --db {params.database} \
        {params.infile} \
        --threads {threads} \
        --output {output.outfmt} \
        --report {output.repfmt} \
        --use-mpa-style \
        --report-zero-counts \
        --use-names  \
        --gzip-compressed
        """

rule kraken_rs:
    output:
        outfmt = "output/{sample}.kraken2.output.rs.tsv",
        repfmt = "output/{sample}.kraken2.report_mpa.rs.tsv"
    message: "Run {wildcards.sample} through kraken2 database w/bacteria, archaea, human"
    params: 
        database = "/mnt/archgen/microbiome_sciences/reference_databases/built/refseq/refseq20191017/kraken2_db/MiniKraken_RefSeq1910/",
        infile = lambda wildcards: SAMPLES[wildcards.sample],
    threads: 32
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate kraken2
        set -u

        kraken2 --db {params.database} \
        {params.infile} \
        --threads {threads} \
        --output {output.outfmt} \
        --report {output.repfmt} \
        --use-mpa-style \
        --report-zero-counts \
        --use-names  \
        --gzip-compressed
        """
