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

reruns = 3
rerun_id=list(range(1,reruns+1))
medioid_id = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
def_radius_id=[0.03, 0.05, 0.07, 0.09]

rule all:
    input: 
            expand("results/{key}/time_combined/time_combined_{rerun_id}_{medioid_id}_{def_radius_id}.tsv", 
                   key=sample_id.keys(), rerun_id=rerun_id, medioid_id=medioid_id, def_radius_id=def_radius_id),
            expand("results/{key}/binbench_results/binbencher_results_{rerun_id}_{medioid_id}_{def_radius_id}.tsv", 
                   key=sample_id.keys(), rerun_id=rerun_id, medioid_id=medioid_id, def_radius_id=def_radius_id),
        # expand("results/{key}/vamb_runs/vamb_from_strobealign_default_params_{rerun_id}_{medioid_id}_{def_radius_id}/vae_clusters_split.tsv", 
        #        key=sample_id.keys(), rerun_id=[1,2,3], medioid_id = [0.05, 0.06, 0.07, 0.08, 0.09, 0.1], def_radius_id=[0.03, 0.05, 0.07, 0.09]) 
# expand("results/{key}/vamb_runs/vamb_from_strobealign_default_params_1/vae_clusters_split.tsv", key=sample_id.keys()),
# params: a = "2"
    # shell:
    #         """
    #         ~/bxc755/miniconda3/bin/parallel --will-cite --dry-run \
    #         '\
    #         echo {{}};
    #         echo hvamb bin default {params.a} --fasta  \
    #         -p  --bamfiles  -m 2000 {{}}; 
    #         '\
    #         ::: a b 
    #         """

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


cores_per_vamb = 10
cores_total = min(127, cores_per_vamb * reruns)
mem_gb_total = min(1990, int((cores_total/128)*1990))
print(f"cores_total:", cores_total, "mem_gb_total:", mem_gb_total)

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
            cores_per_vamb = cores_per_vamb,
            reruns = reruns,
        threads: cores_total
        # log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_{rerun_id}_{medioid_id}_{def_radius_id}" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_total
        conda: "vamb_changed_v5.0.0"
        shell:
            """
            export _MEDOID_RADIUS={wildcards.medioid_id}
            export _DEFAULT_RADIUS={wildcards.def_radius_id}
            rm -rf {params.dir_name}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id}
            vamb bin default --outdir {params.dir_name}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id}  --fasta {input.contig} \
            -p {params.cores_per_vamb} --bamfiles {input.bamfiles} -m 2000 2> {log}_{wildcards.rerun_id}_{wildcards.medioid_id}_{wildcards.def_radius_id} ; 
            """


rulename = "get_time"
rule get_time:
        input: 
            vamb_log = "log/{key}_vamb_for_strobealign_default_{rerun_id}_{medioid_id}_{def_radius_id}", 
        output:
            time = "results/{key}/time_combined/time_combined_{rerun_id}_{medioid_id}_{def_radius_id}.tsv"
        params:
        threads: cores_total
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
        # vamb_bins_minimap = "results/{key}/vamb_from_minimap/vae_clusters_split.tsv",
        # vamb_bins_strobealign_default = "results/{key}/vamb_from_strobealign_default_params/vae_clusters_split.tsv",
        # vamb_bins_aemb = "results/{key}/vamb_from_strobealign_aemb/vae_clusters_split.tsv",
        ref = "/maps/projects/rasmussen/scratch/ptracker/data/refs_circular_ef/{key}.json"
    output: 
        results = "results/{key}/binbench_results/binbencher_results_{rerun_id}_{medioid_id}_{def_radius_id}.tsv"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    shell: 
        """
        cat {input.vamb_bins_minimap} | sed 's/sample_//g' > {input.vamb_bins_minimap}.renamed

        #               Reference, OnlyOrganisms, Bins(which format?), Assembly?
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
        ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
        """

# fast vamb -n 5 -l 2
# MEDIOID 0.06 now
# ::: 0.05 0.06 0.07 0.08 0.09 0.1
# DEFAULT_RADIUS  0.05 
# ::: 0.03 0.05 0.07 0.09 

# ::: $(seq 1 {params.reruns})  
# TODO remove -e 10

            # rm -rf {output.dir}_{{}};

# rulename = "vamb_for_strobealign_default"
# rule vamb_for_strobealign_default:
#         input: 
#             # bamfiles = expand_dir("results/[key]/strobealign_[value].sorted.bam", sample_id),
#             bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
#             contig = "results/{key}/contigs.flt.fna",
#         output:
#             dir = directory("results/{key}/vamb_from_strobealign_default_params"),
#             vamb_bins = "results/{key}/vamb_from_strobealign_default_params/vae_clusters_split.tsv",
#         threads: threads_fn(rulename)
#         log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
#         benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
#         resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#         conda: "vamb_works2"
#         shell:
#             """
#             rm -rf {output.dir} 
#             vamb bin default --outdir {output.dir} --fasta {input.contig} \
#             -p {threads} --bamfiles {input.bamfiles} -m 2000 
#             """
# print(expand("out.{a}", a = [1,2,3]))

# rulename="binbench_minimap"
# rule binbench_minimap:
#     input:
#         vamb_bins_minimap = "results/{key}/vamb_from_minimap/vae_clusters_split.tsv",
#         vamb_bins_strobealign_default = "results/{key}/vamb_from_strobealign_default_params/vae_clusters_split.tsv",
#         vamb_bins_aemb = "results/{key}/vamb_from_strobealign_aemb/vae_clusters_split.tsv",
#         ref = "/maps/projects/rasmussen/scratch/ptracker/data/refs_circular_ef/{key}.json"
#     output: 
#         results = "results/{key}/binbencher_results",
