################################################################################
# Run HUMAnN3 on Pacific calculus and DeepEvo calculus samples
#
# Irina Velsko, 08/10/2021
################################################################################

from glob import glob
import os
import re

workdir: "/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/humann3"

#### SAMPLES ###################################################################
SAMPLES = {}
for sample in glob("/mnt/archgen/microbiome_calculus/pacific_calculus/04-Analysis/humann3/input/*.gz"):
	SAMPLES[os.path.basename(sample).split(".u")[0]] = sample
################################################################################

if not os.path.isdir("snakemake_tmp"):
    os.makedirs(f"{os.getcwd()}/snakemake_tmp")


rule all:
    input: 
        expand("output/{sample}.unmapped_genefamilies.tsv", sample=SAMPLES.keys()),
        "genefamilies_joined.tsv",
        "genefamilies_joined_cpm.tsv",
        "genefamilies_joined_cpm_ur90rxn.tsv",
        "genefamilies_joined_cpm_ur90rxn_names.tsv",
        "genefamilies_joined_cpm_ko.tsv",
        "genefamilies_joined_cpm_ko_names.tsv",
        "pathabundance_joined.tsv",
        "pathabundance_joined_cpm.tsv"       


rule humann2:
    output:
        "output/{sample}.unmapped_genefamilies.tsv"
    message: "Run humann3 on {wildcards.sample}"
    params: 
        fastq = lambda wildcards: SAMPLES[wildcards.sample]
    threads: 24
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann3 --input {params.fastq} --output output --threads {threads}
        """

rule join_gf:
    input:
        "output/"
    output:
        "genefamilies_joined.tsv"
    message: "Run humann3 join on gene families"
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_join_tables -i {input} -o {output} --file_name unmapped_genefamilies
       """

rule renorm_gf:
    input:
        "genefamilies_joined.tsv"
    output:
        "genefamilies_joined_cpm.tsv"
    message: "Run humann3 renorm on gene families"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_renorm_table --input {input} --output {output} --units cpm
       """

rule regroup_gf_ur90:
    input:
        "genefamilies_joined_cpm.tsv"
    output:
        "genefamilies_joined_cpm_ur90rxn.tsv"
    message: "Run humann3 regroup on gene families for UR90"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_regroup_table --input {input} --output {output} --groups uniref90_rxn
       """

rule rename_gf_ur90:
    input:
        "genefamilies_joined_cpm_ur90rxn.tsv"
    output:
        "genefamilies_joined_cpm_ur90rxn_names.tsv"
    message: "Run humann3 rename on gene families for UR90"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_rename_table --input {input} --output {output} -n uniref90
       """


rule regroup_gf_kegg:
    input:
        "genefamilies_joined_cpm.tsv"
    output:
        "genefamilies_joined_cpm_ko.tsv"
    message: "Run humann3 regroup on gene families for KO"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_regroup_table --input {input} --output {output} --groups uniref90_ko
       """

rule rename_gf_ko:
    input:
        "genefamilies_joined_cpm_ko.tsv"
    output:
        "genefamilies_joined_cpm_ko_names.tsv"
    message: "Run humann3 rename on gene families for KO"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_rename_table --input {input} --output {output} -n kegg-orthology
       """

rule join_pa:
    input:
        "output/"
    output:
        "pathabundance_joined.tsv"
    message: "Run humann3 join on path abundance"
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_join_tables -i {input} -o {output} --file_name unmapped_pathabundance
       """

rule renorm_pa:
    input:
        "pathabundance_joined.tsv"
    output:
        "pathabundance_joined_cpm.tsv"
    message: "Run humann3 renorm on gene families"
    params: 
    shell:
        """
        set +u
        source $HOME/miniconda3/etc/profile.d/conda.sh
        conda activate humann3
        set -u
        
        humann_renorm_table --input {input} --output {output} --units cpm
       """

