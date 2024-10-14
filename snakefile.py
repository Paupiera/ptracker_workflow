configfile: "config/config.yaml"

import pandas as pd
import collections
import os

# shell.prefix("source activate ~/bxc755/miniconda3/envs/ptracker_pipeline4; ")
shell.prefix("""
    module unload gcc/13.2.0
    module unload gcc/12.2.0
    module load gcc/13.2.0;
""")

# Paths
OUTDIR= "outdir_plamb" #config["outdir"] #get_config('outdir', 'outdir_plamb', r'.*') # TODO fix
PAU_SRC_DIR = "bin/ptracker/src/workflow"  

# Define deault threads/walltime/mem_gb
default_walltime = config.get("default_walltime")
default_threads = config.get("default_threads")
default_mem_gb = config.get("default_mem_gb")

# Functions to get the config-defined threads/walltime/mem_gb for a rule and if not defined the default
threads_fn = lambda rulename: config.get(rulename, {"threads": default_threads}).get("threads", default_threads) 
walltime_fn  = lambda rulename: config.get(rulename, {"walltime": default_walltime}).get("walltime", default_walltime) 
mem_gb_fn  = lambda rulename: config.get(rulename, {"mem_gb": default_mem_gb}).get("mem_gb", default_mem_gb) 

# Read in the sample data
df = pd.read_csv(config["files"], sep="\s+", comment="#")
sample_id = collections.defaultdict(list)
sample_id_path = collections.defaultdict(dict)
for sample, id, read1, read2 in zip(df.SAMPLE, df.ID, df.READ1, df.READ2):
    id = str(id)
    sample = str(sample)
    sample_id[sample].append(id)
    sample_id_path[sample][id] = [read1, read2]

# Print out run information
print("Running for the following:")
for sample in sample_id.keys():
    print("-"*20)
    print("Sample:", f"{sample}:")
    for id in sample_id[sample]:
        print(f"{id}:")
        print(sample_id_path[sample][id])
    print("-"*20)

#  Define paths to the reads
read_fw  = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][0]
read_rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.id][1]
#  And to the reads after qc
read_fw_after_fastp = "data/sample_{key}/reads_fastp/{id}_1.qc.fastq.gz" 
read_rv_after_fastp =  "data/sample_{key}/reads_fastp/{id}_2.qc.fastq.gz"

## Contig parameters
CONTIGS = config.get("contigs") #get_config('contigs', 'contigs.txt', r'.*') # each line is a contigs path from a given sample
MIN_CONTIG_LEN = int(config.get("min_contig_len")) #get_config('min_contig_len', '2000', r'[1-9]\d*$'))

## N2V parameters
N2V_NZ= config.get("n2v_nz", "weight") 
N2V_ED= config.get("n2v_ed", 128) 
N2V_WL= config.get("n2v_wl", 10) 
N2V_NW= config.get("n2v_nw", 50) 
N2V_WS= config.get("n2v_ws", 10) 
N2V_P= config.get("n2v_p", 0.1) 
N2V_Q= config.get("n2v_q", 2.0) 

NEIGHS_R=config.get("neighs_r", '0.05') 

## Binning parameters
PLAMB_PARAMS = config.get("plamb_params", ' -o C --minfasta 200000  ') 
PLAMB_PRELOAD = config.get("plamb_preload", "") 

# is GPU used ? #
CUDA = config.get("cuda", False)

try: # TODO why?
    os.makedirs(os.path.join(OUTDIR,'log'), exist_ok=True)
except FileExistsError:
    pass

rule all:
    input:
        expand(os.path.join(OUTDIR, "{key}", 'log/run_vamb_asymmetric.finished'), key=sample_id.keys()),
        expand(os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_graph_thr_0.75_candidate_plasmids.tsv'),key=sample_id.keys()),
        expand(os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_unsplit_geNomadplasclustercontigs_extracted_thr_0.75_thrcirc_0.5.tsv'),key=sample_id.keys()),
        expand(os.path.join(OUTDIR,"{key}",'log/run_geNomad.finished'), key=sample_id.keys()),
        # expand("data/sample_{key}/vamb_default", key=sample_id.keys()),
        # expand("data/sample_{key}/vamb_default", key=sample_id.keys()),
        # expand_dir("data/sample_[key]/scapp_[value]/delete_me", sample_id)
        #expand_dir("data/sample_[key]/mp_spades_[value]/contigs.fasta", sample_id),

# rulename = "fastp"
# rule fastp:
#    input: 
#       fw = read_fw, 
#       rv = read_rv, 
#    output:
#       html = "data/sample_{key}/reads_fastp/{id}/report.html", # TODO insert stuff
#       json = "data/sample_{key}/reads_fastp/{id}/report.json",
#       fw = read_fw_after_fastp, 
#       rv = read_rv_after_fastp, 
#    threads: threads_fn(rulename)
#    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
#    log: config.get("log", "log/") + "{key}_{id}_" + rulename
#    shell:
#            'bin/fastp -i {input.fw:q} -I {input.rv:q} '
#            '-o {output.fw:q} -O {output.rv:q} --html {output.html:q} --json {output.json:q} '
#            '--trim_poly_g --poly_g_min_len 7 --cut_tail --cut_front '
#            '--cut_window_size 6  '
#            '--thread {threads} 2> {log:q}'

# rulename = "spades"
# rule spades:
#     input:
#        fw = read_fw_after_fastp, 
#        rv = read_rv_after_fastp, 
#     output:
#        outdir = directory("data/sample_{key}/spades_{id}"),
#        outfile = "data/sample_{key}/spades_{id}/contigs.fasta",
#        graph = "data/sample_{key}/spades_{id}/assembly_graph_after_simplification.gfa", # The graph Changed
#        graphinfo  = "data/sample_{key}/spades_{id}/contigs.paths", # The graph Changed
#     threads: threads_fn(rulename)
#     resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
#     benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
#     log: config.get("log", "log/") + "{key}_{id}_" + rulename
#     shellsmk
#        "bin/SPAdes-3.15.4-Linux/bin/metaspades.py "
#        "-t 20 -m 180 "
#        "-o {output.outdir} -1 {input.fw} -2 {input.rv} " 
#        "-t {threads} --memory {resources.mem_gb} > {log} " 

rulename = "rename_contigs"
rule rename_contigs:
    input:
        "data/sample_{key}/spades_{id}/contigs.fasta"
    output:
        "data/sample_{key}/spades_{id}/contigs.renamed.fasta"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", "log/") + "{key}_{id}_" + rulename
    shell:
        """
        sed 's/^>/>S{wildcards.id}C/' {input} > {output} 2> {log}
        """

rulename="cat_contigs"
rule cat_contigs:
    input: lambda wildcards: expand("data/sample_{key}/spades_{id}/contigs.renamed.fasta", key=wildcards.key, id=sample_id[wildcards.key]),
        # expand_dir("data/sample_[key]/spades_[value]/contigs.renamed.fasta", sample_id)
    output: "data/sample_{key}/contigs.flt.fna.gz"
    threads: threads_fn(rulename)
    params: script = "bin/vamb/src/concatenate.py"
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell: 
        "python {params.script} {output} {input} 2> {log} "  # TODO should filter depending on size????

rulename = "get_contig_names"
rule get_contig_names:
    input:
        "data/sample_{key}/contigs.flt.fna.gz"
    output: 
        "data/sample_{key}/contigs.names.sorted"
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    shell:
        "zcat {input} | grep '>' | sed 's/>//' > {output} 2> {log} "


rulename = "Strobealign_bam_default"
rule Strobealign_bam_default:
        input: 
            fw = read_fw_after_fastp,
            rv = read_rv_after_fastp,
            contig = "data/sample_{key}/contigs.flt.fna.gz",
            # fw = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][0],
            # rv = lambda wildcards: sample_id_path[wildcards.key][wildcards.value][1],
            # contig = "results/{key}/contigs.flt.fna",
        output:
            "data/sample_{key}/mapped/{id}.bam"
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
        log: config.get("log", "log/") + "{key}_{id}_" + rulename
        conda: "envs/strobe_env.yaml"
        shell:
            """
            # module load samtools
            strobealign -t {threads} {input.contig} {input.fw} {input.rv} | samtools sort -o {output} 2> {log}
            """

 # Sort bam files
rulename="sort"
rule sort:
    input:
        "data/sample_{key}/mapped/{id}.bam",
    output:
        "data/sample_{key}/mapped/{id}.bam.sort",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", "log/") + "{key}_{id}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        "samtools sort {input} -o {output} "

## The next part of the pipeline is composed of the following steps:
# (Despite the steps are numerated, some of the order might change)
#   1. Align contigs all against all 
#   2. Generate per sample assembly graphs from gfa assembly graphs
#   3. Run n2v on the per sample assembly graphs 
#   4. Extract hoods from assembly graphs n2v embeddings per sample 
#   5. Connect hoods with alignment graphs
#   6. Run vamb to merge the hoods
#   7. Classify bins/clusters into plasmid/organism/virus bins/clusters

# Align contigs all against all
rulename = "align_contigs"
rule align_contigs:
    input:
        "data/sample_{key}/contigs.flt.fna.gz",
    output:
        os.path.join(OUTDIR,"{key}",'tmp','blastn','blastn_all_against_all.txt'),
        os.path.join(OUTDIR,"{key}",'log/blastn/align_contigs.finished')
    params:
        db_name=os.path.join(OUTDIR,'tmp', "{key}",'blastn','contigs.db'), # TODO should be made?
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    log: config.get("log", "log/") + "{key}_" + rulename
    shell:
        """
        gunzip -c {input} |makeblastdb -in - -dbtype nucl -out {params.db_name} -title contigs.db 
        gunzip -c {input} |blastn -query - -db {params.db_name} -out {output[0]} -outfmt 6 -perc_identity 95 -num_threads {threads} -max_hsps 1000000 >> {log}
        touch {output[1]}
        """
        
## Generate nx graph per sample for graphs from gfa assembly graphs
# for each "sample"
rulename = "weighted_assembly_graphs"
rule weighted_assembly_graphs:
    input:
        graph = "data/sample_{key}/spades_{id}/assembly_graph_after_simplification.gfa", # The graph Changedhis created?
        graphinfo  = "data/sample_{key}/spades_{id}/contigs.paths", # TODO The graph Changed Where is this
    output:
        os.path.join(OUTDIR,"{key}",'tmp','assembly_graphs','{id}.pkl'),
        os.path.join(OUTDIR,"{key}", 'log','assembly_graph_processing','weighted_assembly_graphs_{id}.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'process_gfa.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", "log/") + "{key}_{id}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} --gfa {input[0]} --paths {input[1]} -s {wildcards.id} -m {MIN_CONTIG_LEN}  --out {output[0]} 2> {log}
        touch {output[1]}
        """

# TODO Why does this exist?
rulename = "weighted_assembly_graphs_all_samples"
rule weighted_assembly_graphs_all_samples:
    input: lambda wildcards: expand(os.path.join(OUTDIR,"{key}",'tmp','assembly_graphs','{id}.pkl'),key=wildcards.key, id=sample_id[wildcards.key]),
        # expand_dir(os.path.join(OUTDIR,"[key]",'tmp','assembly_graphs','[value].pkl'), sample_id), 
    output:
        os.path.join(OUTDIR,"{key}",'log','assembly_graph_processing','weighted_assembly_graphs_all_samples.finished')
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    shell:
        """
        touch {output}
        """

## 3. Genereate nx graph from the alignment graph
rulename = "weighted_alignment_graph"
rule weighted_alignment_graph:
    input:
        os.path.join(OUTDIR,"{key}",'tmp','blastn','blastn_all_against_all.txt'),
        os.path.join(OUTDIR,"{key}",'log/blastn/align_contigs.finished')
    output:
        os.path.join(OUTDIR,"{key}",'tmp','alignment_graph','alignment_graph.pkl'),
        os.path.join(OUTDIR,"{key}",'log','alignment_graph_processing','weighted_alignment_graph.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'process_blastout.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} --blastout {input[0]} --out {output[0]} --minid 98 2> {log}
        touch {output[1]}
        """


## 4. Merge assembly graphs an alignment graph into a unified graph
rulename = "create_assembly_alignment_graph"
rule create_assembly_alignment_graph:
    input:
        alignment_graph_file = os.path.join(OUTDIR,"{key}",'tmp','alignment_graph','alignment_graph.pkl'),
        # assembly_graph_files = expand_dir(os.path.join(OUTDIR,"[key]",'tmp','assembly_graphs','[value].pkl'), sample_id), # TODO might be funky 
        assembly_graph_files = lambda wildcards: expand(os.path.join(OUTDIR,"{key}",'tmp','assembly_graphs','{id}.pkl'), key=wildcards.key, id=sample_id[wildcards.key]),
        weighted_alignment_graph_finished_log = os.path.join(OUTDIR,"{key}",'log','alignment_graph_processing','weighted_alignment_graph.finished'),
        weighted_assembly_graphs_all_samples_finished_log = os.path.join(OUTDIR,"{key}", 'log','assembly_graph_processing','weighted_assembly_graphs_all_samples.finished')
    output:
        os.path.join(OUTDIR,"{key}",'tmp','assembly_alignment_graph.pkl'),
        os.path.join(OUTDIR,"{key}", 'log','create_assembly_alignment_graph.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'merge_assembly_alignment_graphs.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} --graph_alignment {input.alignment_graph_file}  --graphs_assembly {input.assembly_graph_files} --out {output[0]}  2> {log}
        touch {output[1]}
        """

## 3. Run n2v on the per sample assembly graphs
rulename = "n2v_assembly_alignment_graph"
rule n2v_assembly_alignment_graph:
    input:
        os.path.join(OUTDIR,"{key}",'tmp','assembly_alignment_graph.pkl'),
        os.path.join(OUTDIR,"{key}",'log','create_assembly_alignment_graph.finished'), 
        contig_names_file = "data/sample_{key}/contigs.names.sorted"
    output:
        directory(os.path.join(OUTDIR,"{key}",'tmp','n2v','assembly_alignment_graph_embeddings')),
        os.path.join(OUTDIR,"{key}",'tmp','n2v','assembly_alignment_graph_embeddings','embeddings.npz'),
        os.path.join(OUTDIR,"{key}",'tmp','n2v','assembly_alignment_graph_embeddings','contigs_embedded.txt'),
        os.path.join(OUTDIR,"{key}",'log','n2v','n2v_assembly_alignment_graph.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'fastnode2vec_args.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} -G {input[0]} --ed {N2V_ED} --nw {N2V_NW} --ws {N2V_WS} --wl {N2V_WL}\
         -p {N2V_P} -q {N2V_Q} --outdirembs {output[0]} --normE {N2V_NZ} --contignames {input.contig_names_file}
        touch {output[3]}
        """


## 4. Extract hoods from assembly graphs n2v embeddings per sample 
rulename = "extract_neighs_from_n2v_embeddings"
rule extract_neighs_from_n2v_embeddings:
    input:
        os.path.join(OUTDIR,"{key}",'tmp','n2v','assembly_alignment_graph_embeddings','embeddings.npz'),
        os.path.join(OUTDIR,"{key}",'tmp','n2v','assembly_alignment_graph_embeddings','contigs_embedded.txt'),
        os.path.join(OUTDIR,"{key}",'log','n2v','n2v_assembly_alignment_graph.finished'),
        os.path.join(OUTDIR,"{key}",'tmp','assembly_alignment_graph.pkl'),
        contig_names_file = "data/sample_{key}/contigs.names.sorted"
    output:
        directory(os.path.join(OUTDIR,"{key}",'tmp','neighs')),
        os.path.join(OUTDIR,"{key}",'tmp','neighs','neighs_object_r_%s.npz'%NEIGHS_R),
        os.path.join(OUTDIR,"{key}",'log','neighs','extract_neighs_from_n2v_embeddings.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'embeddings_to_neighs.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} --embs {input[0]} --contigs_embs {input[1]}\
         --contignames {input.contig_names_file} -g {input[3]} -r {NEIGHS_R} --neighs_outdir {output[0]}
        touch {output[2]}
        """

## 5. Run vamb to merge the hoods
rulename = "run_vamb_asymmetric"
rule run_vamb_asymmetric:
    input:
        notused = os.path.join(OUTDIR,"{key}",'log','neighs','extract_neighs_from_n2v_embeddings.finished'), # TODO why is this not used?
        contigs = "data/sample_{key}/contigs.flt.fna.gz",
        # bamfiles = expand_dir("data/sample_[key]/mapped/[value].bam.sort", sample_id), 
        bamfiles = lambda wildcards: expand("data/sample_{key}/mapped/{id}.bam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
        nb_file = os.path.join(OUTDIR,"{key}",'tmp','neighs','neighs_object_r_%s.npz'%NEIGHS_R)#,
    output:
        directory = directory(os.path.join(OUTDIR,"{key}", 'vamb_asymmetric')),
        bins = os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_unsplit.tsv'),
        finished = os.path.join(OUTDIR,"{key}",'log/run_vamb_asymmetric.finished'),
        lengths = os.path.join(OUTDIR,"{key}",'vamb_asymmetric','lengths.npz'),
    params:
        walltime='86400',
        cuda='--cuda' if CUDA else ''
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        rmdir {output.directory}
        {PLAMB_PRELOAD}
        vamb bin vae_asy --outdir {output.directory} --fasta {input.contigs} -p {threads} --bamfiles {input.bamfiles}\
        --seed 1 --neighs {input.nb_file}  -m {MIN_CONTIG_LEN} {PLAMB_PARAMS}\
         {params.cuda}  
        touch {output}
        """

rulename = "run_geNomad"
rule run_geNomad:
    input:
        #CONTIGS_FILE
        "data/sample_{key}/contigs.flt.fna.gz",
    output:
        directory(os.path.join(OUTDIR,"{key}",'tmp','geNomad')),
        os.path.join(OUTDIR,"{key}",'log/run_geNomad.finished'),
        os.path.join(OUTDIR,"{key}",'tmp','geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv')
    params:
        db_geNomad="genomad_db",
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        genomad end-to-end --cleanup {input} {output[0]}   {params.db_geNomad} --threads {threads}
        touch {output[1]}
        """

## 7. Classify bins/clusters into plasmid/organism/virus bins/clusters
rulename = "classify_bins_with_geNomad"
rule classify_bins_with_geNomad:
    input:
        os.path.join(OUTDIR,"{key}",'tmp','geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv'),
        os.path.join(OUTDIR,"{key}",'log/run_vamb_asymmetric.finished'),
        os.path.join(OUTDIR,"{key}",'log/run_geNomad.finished'),
        os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_unsplit.tsv'),
        contignames = "data/sample_{key}/contigs.names.sorted",
        lengths = os.path.join(OUTDIR,"{key}",'vamb_asymmetric','lengths.npz'),
    output:
        os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_graph_thr_0.75_candidate_plasmids.tsv'),
        os.path.join(OUTDIR,"{key}",'vamb_asymmetric','vae_clusters_unsplit_geNomadplasclustercontigs_extracted_thr_0.75_thrcirc_0.5.tsv'),
        os.path.join(OUTDIR,"{key}",'log','classify_bins_with_geNomad.finished')
    params:
        path = os.path.join(PAU_SRC_DIR, 'src', 'classify_bins_with_geNomad.py'),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    conda: "envs/pipeline_conda.yaml"
    shell:
        """
        python {params.path} --clusters {OUTDIR}/{wildcards.key}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit.tsv\
         --dflt_cls {OUTDIR}/{wildcards.key}/vamb_asymmetric/vae_clusters_unsplit.tsv --scores {input[0]} --outp {output[0]} --lengths {input.lengths} --contignames {input.contignames}
        touch {output[1]}
        """



## EXTRA FOR TESTING 

rulename = "VAMB_DEFAULT"
rule VAMB_DEFAULT:
    input: 
        contig = "data/sample_{key}/contigs.flt.fna.gz",
        bamfiles = lambda wildcards: expand("data/sample_{key}/mapped/{id}.bam.sort", key=wildcards.key, id=sample_id[wildcards.key]),
    output:
        dir = directory("data/sample_{key}/vamb_default"),
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_" + rulename
    log: config.get("log", "log/") + "{key}_" + rulename
    shell:
        """
        rm -rf {output.dir} 
        vamb bin default --outdir {output.dir} --fasta {input.contig} \
        -p {threads} --bamfiles {input.bamfiles} -m 2000 
        """


rulename = "SCAPP"
rule SCAPP:
    input: 
        graph = "data/sample_{key}/spades_{id}/assembly_graph.fastg",
        fw = read_fw_after_fastp, 
        rv = read_rv_after_fastp, 
        # graph_align = "data/sample_{key}/old_scapp/scapp_{id}/assembly_graph.confident_cycs.fasta/intermediate_files/reads_pe_primary.sort.bam", 
    output:
        # scapp makes a directory first and the changes it to a file which is the output. This confuses
        # snakemakes due it either looking for a directory or a file and then stopping the job
        # Therefore need to look for a temp, file which is made when the job is done..
        fake_output = "data/sample_{key}/scapp_{id}/delete_me"
    params:
        true_output = "data/sample_{key}/scapp_{id}/assembly_graph.confident_cycs.fasta", 
    threads: threads_fn(rulename)
    resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
    benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
    log: config.get("log", "log/") + "{key}_{id}_" + rulename
    conda: "envs/install_scapp.yaml"
    shell:
        """
        scapp -g {input.graph} -o {params.true_output} -r1 {input.fw} -r2 {input.rv} -p {threads}
        touch {output.fake_output}
        """
# Expected output file
# fd assembly_graph.confident_cycs.fasta -t f
# scapp_13/assembly_graph.confident_cycs.fasta/assembly_graph.confident_cycs.fasta


rulename = "mpSpades"
rule mpSpades:
        input: 
            fw = read_fw_after_fastp, 
            rv = read_rv_after_fastp, 
        output:
            outdir = directory("data/sample_{key}/mp_spades_{id}"),
            outfile = "data/sample_{key}/mp_spades_{id}/contigs.fasta",
        threads: threads_fn(rulename)
        resources: walltime = walltime_fn(rulename), mem_gb = mem_gb_fn(rulename)
        benchmark: config.get("benchmark", "benchmark/") + "{key}_{id}_" + rulename
        log: config.get("log", "log/") + "{key}_{id}_" + rulename
        shell:
            """
            rm -rf {output.outdir}
            /maps/projects/rasmussen/scratch/ptracker/run_mp_spades/bin/SPades-4/SPAdes-4.0.0-Linux/bin/metaplasmidspades.py --phred-offset 33 -o {output.outdir} -1 {input.fw} -2 {input.rv} \
            > {log}
            ## bin/SPAdes-3.15.4-Linux/bin/metaplasmidspades.py --phred-offset 33 -o {output.outdir} -1 {input.fw} -2 {input.rv} \
            """
