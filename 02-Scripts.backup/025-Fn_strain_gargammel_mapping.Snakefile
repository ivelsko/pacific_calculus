################################################################################
# Map Pacific calculus against Fusobacterium nucleatum reference genomes for polymut,
# to compare the polymut score for mapping to the "wrong" reference strain
#
# Irina Velsko 21/03/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/F_nucleatum_gargammel" 

#### SAMPLES ###################################################################
SAMPLES = [os.path.basename(fn).replace("_s_1.fq.gz", "")
           for fn in glob("/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/sim_shortread_data/art/*_1.fq.gz")
           if os.path.basename(fn).startswith("F_nucleatum")]
#for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/F_nucleatum_gargammel/input/*.fq.gz"):
	#SAMPLES[os.path.basename(sample).split(".")[0]] = sample
################################################################################

#### GENOMES ###################################################################
GENOMES = {'Fnn': 'F_nucleatum_nucleatum',
           'Fnp': 'F_nucleatum_polymorphum',
           'Fnv': 'F_nucleatum_vincentii'}

DEPTH = {'F_nucleatum_nucleatum': (30000000, 2174500),
         'F_nucleatum_polymorphum': (15000000, 2622370),
         'F_nucleatum_vincentii': (15000000, 2268272)}
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam.bai", sample=[s for s in SAMPLES if "short" in s], genome=GENOMES)
        #expand("bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam", sample=SAMPLES, genome=GENOMES),
        #expand("bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam.bai", sample=SAMPLES, genome=GENOMES)
        #expand("bwa_mapping_prokka/{sample}.Fnn.m_s_r.bam", sample=SAMPLES.keys()),
        #expand("bwa_mapping_prokka/{sample}.Fnn.m_s_r.bam.bai", sample=SAMPLES.keys()),
        #expand("bwa_mapping_prokka/{sample}.Fnp.m_s_r.bam", sample=SAMPLES.keys()),
        #expand("bwa_mapping_prokka/{sample}.Fnp.m_s_r.bam.bai", sample=SAMPLES.keys()),
        #expand("bwa_mapping_prokka/{sample}.Fnv.m_s_r.bam", sample=SAMPLES.keys()),
        #expand("bwa_mapping_prokka/{sample}.Fnv.m_s_r.bam.bai", sample=SAMPLES.keys())



rule adapterremoval:
    input:
        sample = ["/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/sim_shortread_data/art/{sample}_s_1.fq.gz",
                  "/mnt/archgen/users/huebner/aDNA-DenovoAssembly/03-data/sim_shortread_data/art/{sample}_s_2.fq.gz"]
    output:
        fq1 = temp("input/{sample}_R1.fastq.gz"),
        fq2 = temp("input/{sample}_R2.fastq.gz"),
        collapsed = temp("input/{sample}.collapsed.fastq.gz"),
        collapsed_trunc = temp("input/{sample}.collapsed_trunc.fastq.gz"),
        singleton = temp("input/{sample}.singleton.fastq.gz"),
        discarded = temp("input/{sample}.discarded.fastq.gz"),
        settings = "input/{sample}.settings" 
    message: "Remove adapters: {wildcards.sample}"
    resources:
        mem = 16,
        cores = 8
    params:
        adapters = "--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGACTATCTCGTATGCCGTCTTCTGCTTG --adapter2 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
        extra = "--collapse --minquality 20 --minlength 30 --trimns --trimqualities --qualitybase 33"
    threads: 8
    wrapper:
        "v1.25.0/bio/adapterremoval"

rule concat:
    input:
        collapsed = "input/{sample}.collapsed.fastq.gz",
        collapsed_trunc = "input/{sample}.collapsed_trunc.fastq.gz",
        singleton = "input/{sample}.singleton.fastq.gz"
    output:
        "input/{sample}.ar.fastq.gz"
    message: "Concatenate FastQ files: {wildcards.sample}"
    resources:
        mem = 2,
        cores = 1
    shell:
        "cat {input} > {output}"

rule subsample:
    input:
        "input/{sample}.ar.fastq.gz"
    output:
        "input/{sample}.seqtk.fastq.gz"
    message: "Subsample FastQ file: {wildcards.sample}"
    resources:
        mem = 8,
        cores = 1
    params:
        n = lambda wildcards: DEPTH[wildcards.sample.split("-")[0]][1] * 3,
        seed = 1001
    wrapper:
        "v1.25.0/bio/seqtk/subsample/se"


rule bwa_sai:
    input:
        "input/{sample}.seqtk.fastq.gz"
    output:
        temp("bwa_mapping_prokka/{sample}/{sample}.{genome}.sai")
    message: "Align {wildcards.sample} against {wildcards.genome} genome using BWA aln"
    params: 
        reffa = lambda wildcards: f"/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/{GENOMES[wildcards.genome]}.fna"
    threads: 4
    shell:
        """
        bwa aln -n 0.01 -l 1024 -t {threads} \
            {params.reffa} \
            {input} > {output}
        """

rule bwa_samse:
    input:
        fq = "input/{sample}.seqtk.fastq.gz",
        sai = "bwa_mapping_prokka/{sample}/{sample}.{genome}.sai"
    output:
        "bwa_mapping_prokka/{sample}.{genome}.m_s.bam"
    message: "Generate alignment file for {wildcards.sample} against {wildcards.genome} genome"
    params:
        reffa = lambda wildcards: f"/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/{GENOMES[wildcards.genome]}.fna"
    shell:
        """
        bwa samse \
            {params.reffa} \
            {input.sai} \
            {input.fq} | \
        samtools view -Sb -F 4 - | samtools sort - -o {output} 
        """

rule samtools_rmdup:
    input:
        "bwa_mapping_prokka/{sample}.{genome}.m_s.bam"
    output:
        "bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam"
    message: "Remove duplicate mapped reads for {wildcards.sample} against {wildcards.genome} genome"
    shell:
        """
        samtools rmdup -s {input} {output}
        """

rule samtools_index:
    input:
        "bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam"
    output:
        "bwa_mapping_prokka/{sample}.{genome}.m_s_r.bam.bai"
    message: "Index bam file for {wildcards.sample} mapped against {wildcards.genome} genome"
    shell:
        """
        samtools index {input}
        """

# map against the F. nucleatum nucleatum genome
#rule bwa_Fnn:
    #output:
        #temp("bwa_mapping_prokka/{sample}/{sample}.Fnn.sai")
    #message: "Align {wildcards.sample} against Fn nucleatum genome using BWA aln"
    #params: 
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_nucleatum.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #threads: 4
    #shell:
        #"""
        #bwa aln -n 0.01 -l 1024 -t {threads} \
            #{params.reffa} \
            #{params.fastq} > {output}
        #"""

#rule bwa_samse_Fnn:
    #input:
        #"bwa_mapping_prokka/{sample}/{sample}.Fnn.sai"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnn.m_s.bam"
    #message: "Generate alignment file for {wildcards.sample} against Fn nucleatum genome"
    #params:
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_nucleatum.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #shell:
        #"""
        #bwa samse \
            #{params.reffa} \
            #{input} \
            #{params.fastq} | \
        #samtools view -Sb -F 4 - | samtools sort - -o {output} 
        #"""

#rule samtools_rmdup_Fnn:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnn.m_s.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnn.m_s_r.bam"
    #message: "Remove duplicate mapped reads for {wildcards.sample} against Fn nucleatum genome"
    #params:
    #shell:
        #"""
        #samtools rmdup -s {input} {output}
        #"""

#rule samtools_index_Fnn:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnn.m_s_r.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnn.m_s_r.bam.bai"
    #message: "Index bam file for {wildcards.sample} mapped against Fn nucleatum genome"
    #params:
    #shell:
        #"""
        #samtools index {input}
        #"""
        
## map against the F. nucleatum polymorphum genome
#rule bwa_Fnp:
    #output:
        #temp("bwa_mapping_prokka/{sample}/{sample}.Fnp.sai")
    #message: "Align {wildcards.sample} against Fn polymorphum reference genome using BWA aln"
    #params: 
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_polymorphum.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #threads: 4
    #shell:
        #"""
        #bwa aln -n 0.01 -l 1024 -t {threads} \
            #{params.reffa} \
            #{params.fastq} > {output}
        #"""

#rule bwa_samse_Fnp:
    #input:
        #"bwa_mapping_prokka/{sample}/{sample}.Fnp.sai"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnp.m_s.bam"
    #message: "Generate alignment file for {wildcards.sample} against Fn polymorphum genome"
    #params:
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_polymorphum.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #shell:
        #"""
        #bwa samse \
            #{params.reffa} \
            #{input} \
            #{params.fastq} | \
        #samtools view -Sb -F 4 - | samtools sort - -o {output} 
        #"""

#rule samtools_rmdup_Fnp:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnp.m_s.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnp.m_s_r.bam"
    #message: "Remove duplicate mapped reads for {wildcards.sample} against Fn polymorphum genome"
    #params:
    #shell:
        #"""
        #samtools rmdup -s {input} {output}
        #"""

#rule samtools_index_Fnp:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnp.m_s_r.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnp.m_s_r.bam.bai"
    #message: "Index bam file for {wildcards.sample} mapped against Fn polymorphum genome"
    #params:
    #shell:
        #"""
        #samtools index {input}
        #"""


## map against the F. nucleatum vincentii genome
#rule bwa_Fnv:
    #output:
        #temp("bwa_mapping_prokka/{sample}/{sample}.Fnv.sai")
    #message: "Align {wildcards.sample} against T. forsythia reference genome using BWA aln"
    #params: 
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_vincentii.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #threads: 4
    #shell:
        #"""
        #bwa aln -n 0.01 -l 1024 -t {threads} \
            #{params.reffa} \
            #{params.fastq} > {output}
        #"""

#rule bwa_samse_Fnv:
    #input:
        #"bwa_mapping_prokka/{sample}/{sample}.Fnv.sai"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnv.m_s.bam"
    #message: "Generate alignment file for {wildcards.sample} against Fn vincentii genome"
    #params:
        #reffa = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_vincentii.fna",
        #fastq = lambda wildcards: SAMPLES[wildcards.sample]
    #shell:
        #"""
        #bwa samse \
            #{params.reffa} \
            #{input} \
            #{params.fastq} | \
        #samtools view -Sb -F 4 - | samtools sort - -o {output} 
        #"""

#rule samtools_rmdup_Fnv:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnv.m_s.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnv.m_s_r.bam"
    #message: "Remove duplicate mapped reads for {wildcards.sample} against Fn vincentii genome"
    #params:
    #shell:
        #"""
        #samtools rmdup -s {input} {output}
        #"""

#rule samtools_index_Fnv:
    #input:
        #"bwa_mapping_prokka/{sample}.Fnv.m_s_r.bam"
    #output:
        #"bwa_mapping_prokka/{sample}.Fnv.m_s_r.bam.bai"
    #message: "Index bam file for {wildcards.sample} mapped against Fn vincentii genome"
    #params:
    #shell:
        #"""
        #samtools index {input}
        #"""
