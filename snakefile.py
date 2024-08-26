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


default_walltime = "01:00:00"
default_threads = 60
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

reruns = 2
rerun_id=list(range(1,reruns+1))
medioid_id = [0.05]
def_radius_id=[0.03]

rule all:
    input: 
            expand("results/{key}/filtered.bam", key=sample_id.keys())
            # expand("results/{key}/msamtools.txt.gz", key=sample_id.keys())


rulename = "coverm"
rule coverm:
        input: 
            fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
            rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
        output:
            "a"
        #     "results/{key}/strobealign_{value}.sorted.bam"
        # threads: threads_fn(rulename)
        # log: return_none_or_default(config, "log", "log/")+"{key}_{value}_" + rulename
        # benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{value}_" + rulename
        # resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        conda: "envs/coverm.yaml"
        shell:
            """
            # coverm
            coverm contig --min-read-aligned-length 80 --min-read-percent-identity \
                95 --min-read-aligned-percent 80
            """

# rulename = "msamtools"
# rule msamtools:
#         input: 
#             bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
#         output:
#             "results/{key}/msamtools.txt.gz",
#         conda: "envs/msamtools.yaml"
#         threads: threads_fn(rulename)
#         resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#         envmodules: "samtools/1.20"
#         shell:
#             """
#                 samtools view {input.bamfiles} \
#                     | msamtools filter -S -bu -l 80 -p 95 -z 80 --besthit - \
#                     | msamtools profile --multi=proportional --label=SAMPLE --unit=ab -o SAMPLE.profile.txt.gz -
#             """

rulename = "msamtools"
rule msamtools:
        input: 
            bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
        output:
            "results/{key}/filtered.bam",
        conda: "envs/msamtools.yaml"
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        shell:
            """
            msamtools filter {input.bamfiles} -bu -l 80 -p 95 -z 80 --besthit > {output}
            """


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


rulename = "vamb_for_strobealign_default"
rule vamb_for_strobealign_default:
        input: 
            bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
            contig = "results/{key}/contigs.flt.fna",
        output:
            vamb_bins = "results/{key}/vamb_runs/vamb_from_strobealign_default_params_{rerun_id}_{medioid_id}_{def_radius_id}/vae_clusters_split.tsv", 
            vamb_log = "log/{key}_vamb_for_strobealign_default_{rerun_id}_{medioid_id}_{def_radius_id}", 
        params: 
            dir_name = directory("results/{key}/vamb_runs/vamb_from_strobealign_default_params"),
            reruns = reruns,
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{rerun_id}_{medioid_id}_{def_radius_id}" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        threads: threads_fn(rulename)
        conda: "vamb_changed_v5.0.0"
        shell:
            """
            export _MEDOID_RADIUS={wildcards.medioid_id}
            export _DEFAULT_RADIUS={wildcards.def_radius_id}
            rm -rf {params.dir_name}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id}
            vamb bin default --outdir {params.dir_name}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id}  --fasta {input.contig} \
            -p {threads} --bamfiles {input.bamfiles} -m 2000 2> {output.vamb_log} ; 
            # -p {threads} --bamfiles {input.bamfiles} -m 2000 2> {log}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id} ; 
            """


rulename = "get_time"
rule get_time:
        input: 
            vamb_log = "log/{key}_vamb_for_strobealign_default_{rerun_id}_{medioid_id}_{def_radius_id}", 
        output:
            time = "results/{key}/time_combined/time_combined_{rerun_id}_{medioid_id}_{def_radius_id}.tsv"
        params:
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        shell:
            """
            # get time for clustering 
            cat {input.vamb_log} | grep -E "Clustered contigs in .* seconds." | cut -d " " -f12,13 \
            | awk '{{print $1, $2, "{wildcards.rerun_id}","{wildcards.medioid_id}", "{wildcards.def_radius_id}"}}' > {output.time}
            """


rulename="binbench_minimap"
rule binbench_minimap:
    input:
        vamb_bins_minimap = "results/{key}/vamb_runs/vamb_from_strobealign_default_params_{rerun_id}_{medioid_id}_{def_radius_id}/vae_clusters_split.tsv", 
        ref = "/maps/projects/rasmussen/scratch/ptracker/data/refs_circular_ef/{key}.json"
    output: 
        results = "results/{key}/binbench_results/binbencher_results_{rerun_id}_{medioid_id}_{def_radius_id}.tsv"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    shell: 
        """
        cat {input.vamb_bins_minimap} | sed 's/sample_//g' > {input.vamb_bins_minimap}.renamed

                      Reference, OnlyOrganisms, Bins(which format?), Assembly?
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        """

rulename="combine_binbench"
rule combine_binbench:
    input:
            binbench = expand("results/{key}/binbench_results/binbencher_results_{rerun_id}_{medioid_id}_{def_radius_id}.tsv", 
                   key=sample_id.keys(), rerun_id=rerun_id, medioid_id=medioid_id, def_radius_id=def_radius_id),
            time = expand("results/{key}/time_combined/time_combined_{rerun_id}_{medioid_id}_{def_radius_id}.tsv",
                   key=sample_id.keys(), rerun_id=rerun_id, medioid_id=medioid_id, def_radius_id=def_radius_id),
    output: 
            "results/binbencher_combined"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    shell: 
            r"""
            ~/bxc755/miniconda3/bin/parallel  --will-cite \
            '\
            cat {{2}} {{2}} {{2}} {{2}}  | \
              paste {{1}} -
            '\
            ::: {input.binbench} :::+ {input.time}
            """

