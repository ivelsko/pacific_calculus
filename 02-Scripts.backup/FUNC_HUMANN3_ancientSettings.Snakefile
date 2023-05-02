####################################################################################################
# Projects: Coprolite evolution
# Part: Composition analysis
# Step: Effect of BowTie2 mapping parameters on HUMAnN3 for ancient samples
#
# HUMAnN3 uses the program BowTie2 twice during its screening analysis in order to infer first
# presence of taxa (via MetaPhlAn3) and then to align the sequencing data against its gene
# catalogue. Using simulation experiments of clipping sequencing data to specific read lengths and
# by processing of paleofeces, we identified that the alignment of short-read sequencing data with
# ancient DNA damage can be improved when removing the requirement of a minimal read length and
# allowing for a single mismatch in the seeds of BowTie2.
#
# In the following, I want to test whether paleofeces also favour from using these ancient DNA
# parameters during the alignment steps that use BowTie2 to recover a richer gene family profile
# that is closer to the ground truth.
#
# modified by Irina Velsko from the script by Alex Huebner, 20/06/22
####################################################################################################

from glob import glob
import gzip
import os
import re

import numpy as np
import pandas as pd

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/humann3"

#### SAMPLES #######################################################################################
# Coprolite samples
SAMPLELIST = pd.read_csv("../../05-Documentation.backup/samples_included.csv",
                         sep="\t", index_col=[0])
CONTSAMPLES = [line.rstrip() for line in open("../../05-Documentation.backup/h3_cont_samples.txt", "rt")]
# SAMPLELIST = SAMPLELIST.loc[~((SAMPLELIST.index.isin(CONTSAMPLES)) | (SAMPLELIST.index.astype(str).str[0] == "N"))]
SAMPLES = {sample: "/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/screening_all_rerun/samtools/filter/{sample}.unmapped.fastq.gz".format(sample=sample)
           for sample in SAMPLELIST.index.values}
# Published paleofeces
# for fn in glob("../03-data/published_paleofeces/*_HG19unmapped.fastq.gz"):
#     sample = os.path.basename(fn).replace("_HG19unmapped.fastq.gz", "")
#     if sample[1:3] == "SM" and sample not in CONTSAMPLES:
#         SAMPLES[sample] = fn
####################################################################################################


rule all:
    input: 
        expand("anc_params_output/{sample}/{sample}_pathcoverage.tsv", sample=SAMPLES.keys())

rule filter_sequences:
    output:
        n = "../tmp/humann3_clippedreads/{sample}.n",
        fastq = "../tmp/humann3_clippedreads/{sample}.fastq.gz"
    message: "Filter reads of sample {wildcards.sample} for minimal length of 50 bp"
    params: 
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    shell:
        """
        bioawk -c fastx '{{
            if (length($seq) >= 50) {{
                print "@" $name;
                print $seq;
                print "+";
                print $qual;
                n += 1
            }}
        }}END{{
            print  n > ("{output.n}")
        }}' {params.fastq} | gzip > {output.fastq}
        """

rule bowtie2:
    input:
        "../tmp/humann3_clippedreads/{sample}.fastq.gz"
    output:
        "../tmp/humann3_clippedreads/{sample}.metaphlan.sam.bz2"
    message: "Align sample {wildcards.sample} against MetaPhlAn3 database"
    params:
        db = "/home/irina_marie_velsko/miniconda3/envs/mpa3/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901"
    threads: 16
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate mpa3
        set -u
        bowtie2 -x {params.db} \
                -U {input} \
                -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --sam-no-hd --sam-no-sq --no-unal \
                --threads {threads} | bzip2 > {output}
        """

rule metaphlan3:
    input:
        sam = "../tmp/humann3_clippedreads/{sample}.metaphlan.sam.bz2",
        nreads = "../tmp/humann3_clippedreads/{sample}.n"
    output: 
        profile = "../04-Analysis/anc_params_output/{sample}.metaphlan.profile.txt"
    message: "Run MetaPhlAn3 with default settings for sample {wildcards.sample}"
    params:
        nreads = lambda wildcards: open(f"../tmp/humann3_clippedreads/{wildcards.sample}.n", "rt").readline().rstrip() if os.path.isfile(f"../tmp/humann3_clippedreads/{wildcards.sample}.n") else 0
    threads: 1
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate mpa3
        set -u
        metaphlan \
            {input.sam} \
            --input_type sam \
            --force \
            --index mpa_v30_CHOCOPhlAn_201901 \
            --ignore_eukaryotes \
            -t rel_ab_w_read_stats \
            --nreads {params.nreads} \
            --sample_id {wildcards.sample} \
            -o {output.profile} \
            --nproc {threads}
        """

rule prepare_metaphlan3_profile:
    input:
        "../04-Analysis/anc_params_output/{sample}.metaphlan.profile.txt"
    output:
        "../tmp/humann3_clippedreads/{sample}.metaphlan.profile.txt"
    message: "Re-format MetaPhlAn3 profile to fit expected input format of HUMAnN3: {wildcards.sample}"
    run:
        with open(output[0], "wt") as outfile:
            for line in open(input[0], "rt"):
                if line.startswith("#mpa_"):
                    outfile.write("#v30\n")
                elif line.startswith("#clade_name"):
                    outfile.write("#clade_name\trelative_abundance\tcoverage\n")
                else:
                    if not line.startswith("#"):
                        fields = line.split("\t")
                        outfile.write(f"{fields[0]}\t{fields[2]}\t{fields[3]}\n")


rule humann3:
    input:
        fastq = "../tmp/humann3_clippedreads/{sample}.fastq.gz",
        profile = "../tmp/humann3_clippedreads/{sample}.metaphlan.profile.txt"
    output:
        "anc_params_output/{sample}/{sample}_pathcoverage.tsv"
    message: "Run HUMAnN3 on sample {wildcards.sample}"
    resources:
        humann3 = 1
    params:
        mp_profile = "../04-Analysis/metaphlan3/{sample}.metaphlan.profile.txt",
        dir = "../04-Analysis/anc_params_output/{sample}"
    log: "../04-Analysis/anc_params_output/{sample}.log"
    threads: 16
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        humann3 -i {input.fastq} \
                -o {params.dir} \
                --o-log {log} \
                --taxonomic-profile {input.profile} \
                --bowtie-options "-D 20 -R 3 -N 1 -L 20 -i S,1,0.50" \
                --remove-temp \
                --threads {threads} \
                --output-basename {wildcards.sample}
        """

################################################################################
