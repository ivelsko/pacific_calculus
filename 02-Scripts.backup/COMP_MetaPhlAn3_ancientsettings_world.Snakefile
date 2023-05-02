################################################################################
# Projects: Pacific dental calculus
# Part: Composition analysis
# Step: Run MetaPhlAn3 with parameters adapted for ancient DNA
#
# MetaPhlAn3 uses the program BowTie2 to align the sequencing data against its
# marker gene catalogue and to infer the compositional profile. By default, two
# parameters are set that are introducing a technical bias against ancient DNA
# sequencing data. First, MetaPhlAn3 uses BowTie2's "very sensitive" mapping
# parameters, which does not allow any mismatches in the seed. This is
# particularly problematic for data with a high levels of ancient DNA damage.
# Second, MetaPhlAn3 by default only uses reads of at least 70 bp.
#
# In the following, I will regenerate the microbial composition profile of the
# dental calculus samples with modified MetaPhlAn3 parameters that allow a
# single mismatch in the seed of BowTie2 and reads as short as 35 bp.
#
# Alex Huebner, 17/06/22
################################################################################

from glob import glob
import os
import pandas as pd


#### SAMPLES ###################################################################
# SAMPLES, = glob_wildcards("/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/screening_all_rerun/samtools/filter/{sample}.unmapped.fastq.gz")
SAMPLES, = glob_wildcards("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/metaphlan3/ancient_parameters/input_bones_blanks/{sample}.unmapped.fastq.gz")
################################################################################

#### Number of reads per sample ################################################
def return_nreads(wildcards):
    nreads = pd.read_csv(checkpoints.summarise_nreads.get(**wildcards).output[0],
                         sep="\t", index_col=[0])
    return nreads.at[wildcards.sample, 'nReads']
################################################################################

localrules: merge_metaphlan3

rule MetaPhlAn3:
    input: 
        expand("tmp/metaphlan3/{sample}.metaphlan.profile.txt", sample=SAMPLES)

rule count_nreads:
    output:
        "tmp/metaphlan3/{sample}.nreads"
    message: "Count the number of reads: {wildcards.sample}"
    params:
        fastq = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/metaphlan3/ancient_parameters/input_bones_blanks/{sample}.unmapped.fastq.gz"
    shell:
        "bioawk -c fastx 'END{{print NR}}' {params.fastq} > {output}"

checkpoint summarise_nreads:
    input:
        expand("tmp/metaphlan3/{sample}.nreads", sample=SAMPLES)
    output:
        "tmp/metaphlan3/summary_nreads.tsv"
    message: "Summarise the number of reads per sample"
    run:
        pd.DataFrame([(os.path.basename(fn).replace(".nreads", ""),
                       int(open(fn, "rt").readline().rstrip()))
                      for fn in input],
                      columns=['sample', 'nReads']) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

rule bowtie2:
    output:
        "tmp/metaphlan3/{sample}.metaphlan.bwt2.sam"
    message: "Align sample {wildcards.sample} against MetaPhlAn3 database"
    params:
        fastq = "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/metaphlan3/ancient_parameters/input_bones_blanks/{sample}.unmapped.fastq.gz",
        db = "/home/irina_marie_velsko/miniconda3/envs/mpa3/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901"
    threads: 16
    shell:
        """
        bowtie2 -x {params.db} \
                -U {params.fastq} \
                -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
                --sam-no-hd --sam-no-sq --no-unal \
                --threads {threads} \
                -S {output} 
        """


rule metaphlan3:
    input:
        bowtie2 = "tmp/metaphlan3/{sample}.metaphlan.bwt2.sam",
        nreads = "tmp/metaphlan3/summary_nreads.tsv"
    output: 
        "tmp/metaphlan3/{sample}.metaphlan.profile.txt"
    message: "Run MetaPhlAn3 with default settings for sample {wildcards.sample}"
    params:
        nreads = lambda wildcards: return_nreads(wildcards)
    threads: 1
    shell:
        """
        metaphlan \
            {input.bowtie2} \
            --input_type sam \
            --force \
            --index mpa_v30_CHOCOPhlAn_201901 \
            --ignore_eukaryotes \
            -t rel_ab_w_read_stats \
            --nreads {params.nreads} \
            --sample_id {wildcards.sample} \
            -o {output} \
            --nproc {threads}
        """
