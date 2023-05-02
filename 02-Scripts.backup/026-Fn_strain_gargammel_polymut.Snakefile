################################################################################
# Polymorphic rate over protein-coding genes with polymut.py for Fusobacterium nucleatum
#
# Irina Velsko 21/03/2023
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/F_nucleatum_gargammel"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/polymut/F_nucleatum_gargammel/bwa_mapping_prokka/*s.bam"):
	SAMPLES[os.path.basename(sample).split(".m")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input:
        expand("polymut_out/{sample}.Fnn.polymut.tsv", sample=SAMPLES.keys()),
        expand("polymut_out/{sample}.Fnp.polymut.tsv", sample=SAMPLES.keys()),
        expand("polymut_out/{sample}.Fnv.polymut.tsv", sample=SAMPLES.keys())
        
rule polymut_Fnn:
    output:
        "polymut_out/{sample}.Fnn.polymut.tsv"
    message: "Run polymut.py on {wildcards.sample} for Fn nucleatum"
    params:
        refgff = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_nucleatum.gff",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate /mnt/archgen/users/huebner/automatic_MAG_refinement/conda/c7f0391e54d836540d3e1f91a1703ac7
        set -u
        
        polymut.py {params.bam} \
        --sortindex \
        --mincov 5 \
        --minqual 30 \
        --dominant_frq_thrsh 0.8 \
        --gff_file {params.refgff} > {output}
        """

rule polymut_Fnp:
    output:
        "polymut_out/{sample}.Fnp.polymut.tsv"
    message: "Run polymut.py on {wildcards.sample} for Fn polymorphum"
    params:
        refgff = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_polymorphum.gff",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate /mnt/archgen/users/huebner/automatic_MAG_refinement/conda/c7f0391e54d836540d3e1f91a1703ac7
        set -u
        
        polymut.py {params.bam} \
        --sortindex \
        --mincov 5 \
        --minqual 30 \
        --dominant_frq_thrsh 0.8 \
        --gff_file {params.refgff} > {output}
        """
rule polymut_Fnv:
    output:
        "polymut_out/{sample}.Fnv.polymut.tsv"
    message: "Run polymut.py on {wildcards.sample} for Fn vincentii"
    params:
        refgff = "/mnt/archgen/microbiome_calculus/pacific_calculus/01-Data/ref_genomes/fusobacterium_nucleatum/prokka/F_nucleatum_vincentii.gff",
        bam = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 8
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate /mnt/archgen/users/huebner/automatic_MAG_refinement/conda/c7f0391e54d836540d3e1f91a1703ac7
        set -u
        
        polymut.py {params.bam} \
        --sortindex \
        --mincov 5 \
        --minqual 30 \
        --dominant_frq_thrsh 0.8 \
        --gff_file {params.refgff} > {output}
        """
