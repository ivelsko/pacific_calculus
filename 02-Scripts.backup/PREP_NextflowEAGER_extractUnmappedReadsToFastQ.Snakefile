################################################################################
# Extract the reads that could not be properly aligned to the human reference
# genome and convert them to FastQ files
#
# Alex Huebner, 23/05/22
# modified 30/05/2022 by Irina Velsko for the Pacific calculus project
################################################################################

from glob import glob
import os

import pandas as pd

#### SAMPLES ###################################################################
samplelist = pd.read_csv("/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/oceania_eager_list_non_collapsed.tsv",
                         sep="\t")
LIBRARIES = samplelist['Library_ID'].unique().tolist()
SAMPLES = samplelist['Sample_Name'].unique().tolist()
################################################################################

#### Auxilliary functions ######################################################

def list_bamfiles(wildcards):
    seqtypes = samplelist.loc[samplelist['Library_ID'] == wildcards.library, 'SeqType'].unique().tolist()
    return " ".join([f"/mnt/archgen/microbiome_calculus/pacific_calculus/03-Preprocessing/eager2_non_collapsed/mapping/bwa/{wildcards.library}_{seqtype}.mapped.bam"
                     for seqtype in seqtypes])

################################################################################

rule all:
    input:
        "05-results.backup/PREP_Nextflow_EAGER_noReads.tsv",
        expand("03-Preprocessing/eager2_non_collapsed/fastqs/per_sample/{sample}.concatenated", sample=SAMPLES),
        expand("03-Preprocessing/eager2_non_collapsed/fastqs/collapsed/{sample}.collapsed", sample=SAMPLES)


#### Extract and merge non-aligned reads #######################################

rule concat_seqtypes:
    output:
        pipe("tmp/eager_extract_unmapped/{library}.concat.bam")
    message: "Concatenate all BAM files: {wildcards.library}"
    conda: "ENVS_samtools.yaml"
    params:
        bam = lambda wildcards: list_bamfiles(wildcards)
    shell:
        "samtools cat -o {output} {params.bam}"

rule samtools_sort_by_name:
    input:
        "tmp/eager_extract_unmapped/{library}.concat.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{library}.nsorted.bam")
    message: "Sort the BAM file by name: {wildcards.library}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 12,
        cores = 4
    threads: 4
    shell:
        """
        samtools sort -T /tmp/{wildcards.library} -@ {threads} -n -o {output} {input}
        """

rule samtools_fixmate:
    input:
        "tmp/eager_extract_unmapped/{library}.nsorted.bam"
    output:
        pipe("tmp/eager_extract_unmapped/{library}.fixmate.bam")
    message: "Fix mate flags: {wildcards.library}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 4
    threads: 1
    shell:
        """
        samtools fixmate -mcu {input} {output}
        """

rule extract_unmapped_reads:
    input:
        "tmp/eager_extract_unmapped/{library}.fixmate.bam"
    output:
        pe1 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_1.fastq.gz",
        pe2 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_2.fastq.gz",
        pe0 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_0.fastq.gz",
    message: "Extract all reads for which are not aligned in a proper pair and convert to fastq: {wildcards.library}"
    conda: "ENVS_samtools.yaml"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -uh -e '(flag.paired && (flag.unmap || flag.munmap)) || (!flag.paired && flag.unmap)' {input} | \
        samtools fastq -1 {output.pe1} \
                       -2 {output.pe2} \
                       -0 {output.pe0} -
        """

rule count_reads:
    input:
        pe1 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_1.fastq.gz",
        pe2 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_2.fastq.gz",
        pe0 = "03-Preprocessing/eager2_non_collapsed/fastqs/{library}_0.fastq.gz"
    output:
        temp("03-Preprocessing/eager2_non_collapsed/fastqs/{library}.n")
    message: "Count the number of reads: {wildcards.library}"
    conda: "ENVS_bioawk.yaml"
    resources:
        mem = 2
    shell:
        """
        if [[ $(gzip -l {input.pe1} | awk 'NR == 2 {{print $1}}') -le 29 ]]; then
            reads_PE1=0
            reads_PE2=0
        else
            reads_PE1=$(bioawk -c fastx 'END{{print NR}}' {input.pe1})
            reads_PE2=$(bioawk -c fastx 'END{{print NR}}' {input.pe2})
        fi
        if [[ $(gzip -l {input.pe0} | awk 'NR == 2 {{print $1}}') -le 29 ]]; then
            reads_PE0=0
        else
            reads_PE0=$(bioawk -c fastx 'END{{print NR}}' {input.pe0})
        fi
        echo -e "{wildcards.library}\t${{reads_PE1}}\t${{reads_PE2}}\t${{reads_PE0}}" > {output}
        """

rule summarise_count_reads:
    input:
        expand("03-Preprocessing/eager2_non_collapsed/fastqs/{library}.n", library=LIBRARIES)
    output:
        "05-results.backup/PREP_Nextflow_EAGER_noReads.tsv"
    message: "Summarise the number of reads per library"
    run:
        pd.concat([pd.read_csv(fn, sep="\t", header=None, names=['sample', 'R1', 'R2'])
                   for fn in input]) \
            .sort_values(['sample']) \
            .to_csv(output[0], sep="\t", index=False)

rule concatenate_samples:
    input:
        lambda wildcards: [f"03-Preprocessing/eager2_non_collapsed/fastqs/{library}.n" for library in samplelist.loc[samplelist['Sample_Name'] == wildcards.sample]['Library_ID'].unique().tolist()]
    output:
        touch("03-Preprocessing/eager2_non_collapsed/fastqs/per_sample/{sample}.concatenated")
    message: "Concatenate the FastQ files across all libraries: {wildcards.sample}"
    params:
        dir = "03-Preprocessing/eager2_non_collapsed/fastqs/per_sample"
    run:
        libs_sqt = samplelist.loc[samplelist['Sample_Name'] == wildcards.sample][['Library_ID', 'SeqType']] \
            .drop_duplicates()
        if "PE" in libs_sqt['SeqType'].tolist():
            libs = libs_sqt.loc[libs_sqt['SeqType'] == "PE"]['Library_ID'].tolist()
            for i in [1, 2]:
                with open(f"{params.dir}/{wildcards.sample}_{i}.fastq.gz", "wb") as outfile:
                    for lib in libs:
                        for line in open(f"{os.path.dirname(params.dir)}/{lib}_{i}.fastq.gz", "rb"):
                            outfile.write(line)
        if "SE" in libs_sqt['SeqType'].tolist():
            libs = libs_sqt.loc[libs_sqt['SeqType'] == "SE"]['Library_ID'].tolist()
            with open(f"{params.dir}/{wildcards.sample}_0.fastq.gz", "wb") as outfile:
                for lib in libs:
                    for line in open(f"{os.path.dirname(params.dir)}/{lib}_0.fastq.gz", "rb"):
                            outfile.write(line)

rule collapse:
    input:
        "03-Preprocessing/eager2_non_collapsed/fastqs/per_sample/{sample}.concatenated"
    output:
        touch("03-Preprocessing/eager2_non_collapsed/fastqs/collapsed/{sample}.collapsed")
    message: "Merge overlapping reads using fastp: {wildcards.sample}"
    conda: "ENVS_fastp.yaml"
    params:
        dir = "03-Preprocessing/eager2_non_collapsed/fastqs/per_sample",
        outdir = "03-Preprocessing/eager2_non_collapsed/fastqs/collapsed"
    shell:
        """
        if [[ -f "{params.dir}/{wildcards.sample}_1.fastq.gz" ]]; then
            fastp --in1 {params.dir}/{wildcards.sample}_1.fastq.gz \
                  --in2 {params.dir}/{wildcards.sample}_2.fastq.gz \
                  --merge \
                  --out1 {params.outdir}/{wildcards.sample}_1.fastq.gz \
                  --out2 {params.outdir}/{wildcards.sample}_2.fastq.gz \
                  --merged_out {params.outdir}/{wildcards.sample}_collapsed.fastq.gz \
                  --overlap_len_require 11 \
                  -A -G -Q -j {params.outdir}/{wildcards.sample}_PE.fastp.json -h /dev/null
        else
            touch {output}
        fi
        """

################################################################################
