configfile: "config/config.yaml"

from utils import expand_dir, populate_dict_of_lists, config_dict
import pandas as pd
import re
import os
import sys
import numpy as np
import collections
# lambda wildcards: expand("results/{key}/{value}_renamed_contig.fasta", key=wildcards.key, value=sample_id[wildcards.key])

# shell.prefix(f"source activate ~/bxc755/miniconda3/envs/strobealign; ")
shell.prefix(f"module unload gcc/13.2.0; module unload gcc/12.2.0; module load gcc/13.2.0;")

def return_none_or_default(dic, key, default_value):
    if dic == None:
            return default_value
    if dic.get(key) == None:
            return default_value
    return dic[key]


default_walltime = "48:00:00"
default_threads = 10
default_mem_gb = 100 

threads_fn = lambda rulename: return_none_or_default(config.get(rulename), "threads", default_threads) 
walltime_fn = lambda rulename: return_none_or_default(config.get(rulename), "walltime", default_walltime) 
mem_gb_fn = lambda rulename: return_none_or_default(config.get(rulename), "mem_gb", default_mem_gb) 

# Read in the sample data
df = pd.read_csv(config["files"], sep="\s+", comment="#")
sample_id = {}
sample_id_path = collections.defaultdict(dict)
for sample, id, read1, read2, contig in zip(df.SAMPLE, df.ID, df.READ1, df.READ2, df.CONTIG):
    id = str(id)
    sample = str(sample)
    populate_dict_of_lists(sample_id, sample, id)
    sample_id_path[sample][id] = [read1, read2, contig]

# Print out run information
print("Running for the following:")
for sample in sample_id.keys():
    print("-"*20)
    print("Sample:", f"{sample}:")
    for id in sample_id[sample]:
        print(f"{id}:") 
        print(sample_id_path[sample][id])
    print("-"*20)

# def abc(string, sample_id, wildcards):

rule all:
    input: 
        # expand_dir("results/[key]/[value]_headers.tsv", sample_id),
        # expand_dir("data/sample_[key]/mapped/minimap_[value].bam", sample_id),
        # expand_dir("results/[key]/strobealign_[value].sorted.bam", sample_id),
        # "results/Airways/vamb_out",
        # expand_dir("results/[key]/abundances_[value].paf", sample_id),
        expand("results/{key}/vamb_from_minimap/vae_clusters_split.tsv", key=sample_id.keys()),
        expand("results/{key}/vamb_from_strobealign_aemb", key=sample_id.keys()),
        expand("results/{key}/vamb_from_strobealign_default_params", key=sample_id.keys()),
        #expand("results/{key}/abundances_collected.npz", key=sample_id.keys()),
        expand("results/{key}/binbencher_results", key=sample_id.keys()),


rulename = "Rename_Contigs"
rule Rename_Contigs:
        input: 
            contig = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][2],
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        output:
            "results/{key}/{value}_renamed_contig.fasta",
        shell:
            """
            cat {input.contig} | sed "/^>/ s/^>NODE/>S{wildcards.value}CNODE/" > {output}
            """

rulename="cat_contigs"
rule cat_contigs:
    input:
        lambda wildcards: expand("results/{key}/{value}_renamed_contig.fasta", key=wildcards.key, value=sample_id[wildcards.key])
        # expand_dir("results/[key]/[value]_renamed_contig.fasta", sample_id)
    output:
        "results/{key}/contigs.flt.fna"
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    conda: "vamb_works2"
    shell: 
        """
        # Needs to use full path for python see snakemake bug: https://github.com/snakemake/snakemake/issues/2861
        ~/bxc755/miniconda3/envs/vamb_works2/bin/python bin/vamb/src/concatenate.py {output} {input} --keepnames -m 2000 --nozip
        """

rulename = "Strobealign_x"
rule Strobealign_x:
        input: 
            fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
            rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
            contig = "results/{key}/contigs.flt.fna",
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
        conda: "strobe_env.yaml"
        output:
            "results/{key}/abundances_{value}.paf"
        shell:
            """
            module load samtools
            strobealign -x -t {threads} {input.contig} {input.fw} {input.rv} > {output} 2> {log} 
            """

rulename = "Strobealign"
rule Strobealign:
        input: 
            fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
            rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
            contig = "results/{key}/contigs.flt.fna",
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
        conda: "strobe_env.yaml"
        output:
            "results/{key}/abundances_{value}.tsv"
        shell:
            """
            module load samtools
            strobealign -t {threads} --aemb {input.contig} {input.fw} {input.rv} > {output} 2> {log}
            """

rulename = "Stobealign_output_to_npz"
rule Strobealing_output_to_npz:
        input: 
            # expand_dir("results/[key]/abundances_[value].tsv", sample_id),
            lambda wildcards: expand("results/{key}/abundances_{value}.tsv", key=wildcards.key, value=sample_id[wildcards.key])
        output:
            npz = "results/{key}/abundances_collected.npz",
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        conda: "strobe_env.yaml"
        shell:
            """
            python bin/create_numpy_abundance.py {input} {output.npz}
            """
# bash bin/format_abundance.sh {input} > {output}.names # is this even used by the new vamb input? 
    
rulename = "Strobealign_bam_default"
rule Strobealign_bam_default:
        input: 
            fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
            rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
            contig = "results/{key}/contigs.flt.fna",
        output:
            "results/{key}/strobealign_{value}.sorted.bam"
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        conda: "strobe_env.yaml"
        shell:
            """
            module load samtools
            strobealign -t {threads} {input.contig} {input.fw} {input.rv} | samtools sort -o {output} 2> {log}
            """


#### MINIMAP ####
# Index resulting contig-file with minimap2
rulename="index"
rule index:
    input:
        contigs = "results/{key}/contigs.flt.fna",
    output:
        mmi = "results/{key}/contigs.flt.mmi"
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    conda: "ptracker_pipeline4"
    shell:
        "minimap2 -d {output} {input} 2> {log} "

# This rule creates a SAM header from a FASTA file.
# We need it because minimap2 for truly unknowable reasons will write
# SAM headers INTERSPERSED in the output SAM file, making it unparseable.
# To work around this mind-boggling bug, we remove all header lines from
# minimap2's SAM output by grepping, then re-add the header created in this rule.
rulename="dict"
rule dict:
    input:
        contigs = "results/{key}/contigs.flt.fna",
    output:
        "results/{key}/contigs.flt.dict"
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    conda: "ptracker_pipeline4"
    shell:
        "module load samtools; "
        "samtools dict {input.contigs} 2> {log} | cut -f1-3 > {output}"

# Generate bam files 
LONG_READS = False
rulename="minimap"
rule minimap:
    input:
        fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
        rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
        mmi = "results/{key}/contigs.flt.mmi",
        dict = "results/{key}/contigs.flt.dict",
    output:
        bam = "results/{key}/minimap_{value}.bam",
        bam_sorted = "results/{key}/minimap_{value}.sorted.bam"
    params:
        long_or_short_read = 'map-pb -L' if LONG_READS else 'sr',
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    conda: "ptracker_pipeline4"
    shell:
        # See comment over rule "dict" to understand what happens here
        "module load samtools; "
        "minimap2 -t {threads} -ax {params.long_or_short_read} {input.mmi} {input.fw} {input.rv} "
        " | grep -v '^@'"
        " | cat {input.dict} - "
        " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
        " > {output.bam}; "
        "samtools sort {output.bam} -o {output.bam_sorted} "


rulename="vamb_npz"
rule vamb_npz:
    input:
        npz = "results/{key}/abundances_collected.npz",
        contig = "results/{key}/contigs.flt.fna"
    output:
        dir = directory("results/{key}/vamb_from_strobealign_aemb"),
        vamb_bins_aemb = "results/{key}/vamb_from_strobealign_aemb/vae_clusters_split.tsv",
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    conda: "vamb_works2"
    shell: 
        """
        rm -rf {output.dir} 
        vamb bin default --outdir {output.dir} --fasta {input.contig} \
        -p {threads} --rpkm {input.npz} -m 2000 --norefcheck
        """

rulename="vamb_minimap"
rule vamb_minimap:
    input:
        bamfiles = lambda wildcards: expand("results/{key}/minimap_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
        # bamfiles = expand_dir("results/[key]/minimap_[value].sorted.bam", sample_id),
        contig = "results/{key}/contigs.flt.fna",
    output:
        vamb_bins = "results/{key}/vamb_from_minimap/vae_clusters_split.tsv",
    threads: threads_fn(rulename)
    params: dir = directory("results/{key}/vamb_from_minimap"),
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    conda: "vamb_works2"
    shell: 
        """
        rm -rf {params.dir} 
        vamb bin default --outdir {params.dir} --fasta {input.contig} \
        -p {threads} --bamfiles {input.bamfiles} -m 2000 
        """

rulename="binbench_minimap"
rule binbench_minimap:
    input:
        vamb_bins_minimap = "results/{key}/vamb_from_minimap/vae_clusters_split.tsv",
        vamb_bins_strobealign_default = "results/{key}/vamb_from_strobealign_default_params/vae_clusters_split.tsv",
        vamb_bins_aemb = "results/{key}/vamb_from_strobealign_aemb/vae_clusters_split.tsv",
        ref = "/maps/projects/rasmussen/scratch/ptracker/data/refs_circular_ef/{key}.json"
    output: 
        results = "results/{key}/binbencher_results",
        # bins_renamed_minimap = "results/{key}/vamb_from_minimap/bins_renamed.txt",
        # bins_renamed_strob_def = "results/{key}/vamb_from_strobealign_default_params/bins_renamed.txt",
        # bins_renamed_aemb = "results/{key}/vamb_from_strobealign_aemb/bins_renamed.txt",
    threads: threads_fn(rulename)
    log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
    shell: 
        """
        cat {input.vamb_bins_minimap} | sed 's/sample_//g' > {input.vamb_bins_minimap}.renamed
        cat {input.vamb_bins_strobealign_default} | sed 's/sample_//g' > {input.vamb_bins_strobealign_default}.renamed
        cat {input.vamb_bins_aemb} | sed 's/sample_//g' > {input.vamb_bins_aemb}.renamed

        #               Reference, OnlyOrganisms, Bins(which format?), Assembly?

        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_strobealign_default}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_aemb}.renamed    true          >> {output.results}

        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_strobealign_default}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_aemb}.renamed    true          >> {output.results}

        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_strobealign_default}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_aemb}.renamed    false          >> {output.results}

        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_strobealign_default}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_aemb}.renamed    false          >> {output.results}
        """

rulename = "vamb_for_strobealign_default"
rule vamb_for_strobealign_default:
        input: 
            # bamfiles = expand_dir("results/[key]/strobealign_[value].sorted.bam", sample_id),
            bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
            contig = "results/{key}/contigs.flt.fna",
        output:
            dir = directory("results/{key}/vamb_from_strobealign_default_params"),
            vamb_bins = "results/{key}/vamb_from_strobealign_default_params/vae_clusters_split.tsv",
        threads: threads_fn(rulename)
        log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        conda: "vamb_works2"
        shell:
            """
            rm -rf {output.dir} 
            vamb bin default --outdir {output.dir} --fasta {input.contig} \
            -p {threads} --bamfiles {input.bamfiles} -m 2000 
            """
