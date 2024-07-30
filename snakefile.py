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

rule all:
    input: 
        expand("results/{key}/vamb_from_strobealign_default_params_1/vae_clusters_split.tsv", key=sample_id.keys()),
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


reruns = 5
cores_per_vamb = 5
rulename = "vamb_for_strobealign_default"
cores_total = min(128, cores_per_vamb * reruns)
mem_gb_total = min(1990, int((cores_total/128)*1990))
print(f"cores_total:", cores_total, "mem_gb_total:", mem_gb_total)

rule vamb_for_strobealign_default:
        input: 
            # bamfiles = expand_dir("results/[key]/strobealign_[value].sorted.bam", sample_id),
            bamfiles = lambda wildcards: expand("results/{key}/strobealign_{value}.sorted.bam", key=wildcards.key, value=sample_id[wildcards.key]),
            contig = "results/{key}/contigs.flt.fna",
        output:
            vamb_bins = expand("results/{key}/vamb_from_strobealign_default_params_{run_id}/vae_clusters_split.tsv", run_id=list(range(1,reruns+1))),
        params: 
            dir_name = directory("results/{key}/vamb_from_strobealign_default_params"),
            cores_per_vamb = cores_per_vamb,
            reruns = reruns,
        threads: cores_total
        log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
        benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_total
        conda: "vamb_works2"
        shell:
            """
            ~/bxc755/miniconda3/bin/parallel --will-cite \
            '\
            vamb bin default --outdir {params.dir_name}_{{}} --fasta {input.contig} \
            -p {params.cores_per_vamb} --bamfiles {input.bamfiles} -m 2000 2> {log}{{}} ; 
            '\
            ::: $(seq 1 {params.reruns}) 
            """

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
#         # bins_renamed_minimap = "results/{key}/vamb_from_minimap/bins_renamed.txt",
#         # bins_renamed_strob_def = "results/{key}/vamb_from_strobealign_default_params/bins_renamed.txt",
#         # bins_renamed_aemb = "results/{key}/vamb_from_strobealign_aemb/bins_renamed.txt",
#     threads: threads_fn(rulename)
#     log: return_none_or_default(config, "log", "log/")+"{key}_" + rulename
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: return_none_or_default(config, "benchmark", "benchmark/")+"{key}_" + rulename
#     shell: 
#         """
#         cat {input.vamb_bins_minimap} | sed 's/sample_//g' > {input.vamb_bins_minimap}.renamed
#         cat {input.vamb_bins_strobealign_default} | sed 's/sample_//g' > {input.vamb_bins_strobealign_default}.renamed
#         cat {input.vamb_bins_aemb} | sed 's/sample_//g' > {input.vamb_bins_aemb}.renamed
#
#         #               Reference, OnlyOrganisms, Bins(which format?), Assembly?
#
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_strobealign_default}.renamed    true          >> {output.results}
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_aemb}.renamed    true          >> {output.results}
#
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    true          >> {output.results}
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_strobealign_default}.renamed    true          >> {output.results}
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_aemb}.renamed    true          >> {output.results}
#
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_strobealign_default}.renamed    false          >> {output.results}
#         ./Binbench.jl {input.ref}  true          {input.vamb_bins_aemb}.renamed    false          >> {output.results}
#
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_minimap}.renamed    false          >> {output.results}
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_strobealign_default}.renamed    false          >> {output.results}
#         ./Binbench.jl {input.ref}  false          {input.vamb_bins_aemb}.renamed    false          >> {output.results}
#         """
