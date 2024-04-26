import re
import os
import sys
import numpy as np
SNAKEDIR = os.path.dirname(workflow.snakefile)
sys.path.append(os.path.join(SNAKEDIR, 'src'))


## This snakemake workflow executes the second part of the  p-tracker pipeline.
# Inputs:
# - Contig fasta file with all contigs from all samples
# - Assembly graphs in GFA format
# - geNomad contig scores
# Outputs (final):
# - Plasmid bins
# - Organism bins
# - Virus bins

## The pipeline is composed of the following steps:
# 1. Align contigs all against all 
# 2. Generate per sample assembly graphs from gfa assembly graphs
# 3. Run n2v on the per sample assembly graphs 
# 4. Extract hoods from assembly graphs n2v embeddings per sample 
# 5. Connect hoods with alignment graphs
# 6. Run vamb to merge the hoods
# 7. Classify bins/clusters into plasmid/organism/virus bins/clusters

# Despite the steps are numerated, some of the order might change 

def get_config(name, default, regex):
    res = config.get(name, default).strip()
    m = re.match(regex, res)
    if m is None:
        raise ValueError(
            f'Config option \'{name}\' is \'{res}\', but must conform to regex \'{regex}\'')
    return res


# set configurations
## Contig parameters
CONTIGS_FILE=get_config('contigs_filtered_file','', r'.*')
ABUNDANCE_FILE=get_config('abundance_file','', r'.*')
CONTIGNAMES_FILE = get_config('contignames_file', '', r'.*') # file containing the sorted contignames
CONTIGS = get_config('contigs', 'contigs.txt', r'.*') # each line is a contigs path from a given sample
MIN_CONTIG_LEN = int(get_config('min_contig_len', '2000', r'[1-9]\d*$'))
#MIN_IDENTITY = float(get_config('min_identity', '0.95', r'.*'))
#SAMPLE_DATA 

## Assembly graph parameters
ASSEMBLY_GRAPHS_DIR=get_config('assembly_graph_dir','', r'.*')

## N2V parameters
N2V_NZ=get_config('n2v_nz','weight', r'.*')
N2V_ED=get_config('n2v_ed','128', r'.*')
N2V_WL=get_config('n2v_wl','10', r'.*')
N2V_NW=get_config('n2v_nw','50', r'.*')
N2V_WS=get_config('n2v_ws','10', r'.*')
N2V_P=get_config('n2v_p','0.1', r'.*')
N2V_Q=get_config('n2v_q','2.0', r'.*')

NEIGHS_R='0.05'
## Binning parameters
PLAMB_MEM = get_config('plamb_mem', '20gb', r'[1-9]\d*GB$')
PLAMB_PPN = get_config('plamb_ppn', '10', r'[1-9]\d*(:gpus=[1-9]\d*)?$')
PLAMB_PARAMS = get_config('plamb_params',' -o C --minfasta 200000  ', r'.*')
PLAMB_PRELOAD = get_config('plamb_preload', '', r'.*')

RENAMING_FILE = get_config('renaming_file', '', r'.*')
OUTDIR= get_config('outdir', 'outdir_plamb', r'.*')

try:
    os.makedirs(os.path.join(OUTDIR,'log'), exist_ok=True)
except FileExistsError:
    pass


# parse if GPUs is needed #
plamb_threads, sep, plamb_gpus = PLAMB_PPN.partition(':gpus=')
PLAMB_PPN = plamb_threads
CUDA = len(plamb_gpus) > 0

## read in sample information ##
SAMPLES=get_config('samples',"",r'.*')

# target rule
rule all:
    input:
        #os.path.join(OUTDIR,'log/run_vamb_asymmetric.finished')
        os.path.join(OUTDIR,'log/rename_clusters_and_hoods.finished')
        #os.path.join(OUTDIR,'log','classify_bins_with_geNomad.finished')

#gunzip -c {input} |blastn -query - -db {output[0]} -out {output[1]} -outfmt 6 -perc_identity 95 -num_threads {threads} -max_hsps 1000000
# 1. Align contigs all against all 
rule align_contigs:
    input:
        CONTIGS_FILE # Just a contig file
    output:
        os.path.join(OUTDIR,'tmp','blastn','blastn_all_against_all.txt'),
        os.path.join(OUTDIR,'log/blastn/align_contigs.finished')
    params:
        db_name=os.path.join(OUTDIR,'tmp','blastn','contigs.db'),
        walltime='86400',
        nodes='1',
        ppn='8'
    resources:
        mem='15GB'
    threads:
        8
    log:
        o = os.path.join(OUTDIR,'qsub/blastn/align_contigs.o'),
        e = os.path.join(OUTDIR,'qsub/blastn/align_contigs.e')
    conda:
        #'envs/blast.yaml'
        'vamb_n2v_master'
    shell:
        """
        module load ncbi-blast/2.15.0+
        gunzip -c {input} |makeblastdb -in - -dbtype nucl -out {params.db_name} -title contigs.db
        gunzip -c {input} |blastn -query - -db {params.db_name} -out {output[0]} -outfmt 6 -perc_identity 95 -num_threads {threads} -max_hsps 1000000
        touch {output[1]}
        """
        
## 2. Generate nx graph per sample for graphs from gfa assembly graphs
# for each "sample"
rule weighted_assembly_graphs:
    input:
        os.path.join(ASSEMBLY_GRAPHS_DIR,'{sample}','assembly_graph_after_simplification.gfa'), # The graph
        os.path.join(ASSEMBLY_GRAPHS_DIR,'{sample}','contigs.paths'), # information from the graph
    output:
        os.path.join(OUTDIR,'tmp','assembly_graphs','{sample}.pkl'),
        os.path.join(OUTDIR,'log','assembly_graph_processing','weighted_assembly_graphs_{sample}.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'process_gfa.py'),
        walltime='86400',
        nodes='1',
        ppn='4'
    resources:
        mem = '4GB'
    threads:
        4
    log:
        o=os.path.join(OUTDIR,'qsub','assembly_graph_processing','weighted_assembly_graphs_{sample}.out'),
        e=os.path.join(OUTDIR,'qsub','assembly_graph_processing','weighted_assembly_graphs_{sample}.err'),
    conda:
        'vamb_n2v_master'
        #'envs/vamb.yaml'
    shell:
        """
        python {params.path} --gfa {input[0]} --paths {input[1]} -s {wildcards.sample} -m {MIN_CONTIG_LEN}  --out {output[0]}
        touch {output[1]}
        """

def weighted_assembly_graphs_all_samples_f(wildcards):
    samples = [sample for sample in SAMPLES.split()]
    weighted_assembly_graphs_all_samples_paths=expand(
        os.path.join(OUTDIR,'tmp','assembly_graphs','{sample}.pkl'),
        sample=samples)
    return weighted_assembly_graphs_all_samples_paths

rule weighted_assembly_graphs_all_samples:
    input:
        weighted_assembly_graphs_all_samples_f
    output:
        os.path.join(OUTDIR,'log','assembly_graph_processing','weighted_assembly_graphs_all_samples.finished')
    params:
        walltime='86400',
        nodes='1',
        ppn='1'
    resources:
        mem = '1GB'
    threads:
        1
    log:
        o=os.path.join(OUTDIR,'qsub','assembly_graph_processing','weighted_assembly_graphs_all_samples.out'),
        e=os.path.join(OUTDIR,'qsub','assembly_graph_processing','weighted_assembly_graphs_all_samples.err')
    conda:
        'vamb_n2v_master'
        #'envs/vamb.yaml'
    shell:
        """
        touch {output}
        """



## 3. Genereate nx graph from the alignment graph
rule weighted_alignment_graph:
    input:
        os.path.join(OUTDIR,'tmp','blastn','blastn_all_against_all.txt'),
        os.path.join(OUTDIR,'log/blastn/align_contigs.finished')
    output:
        os.path.join(OUTDIR,'tmp','alignment_graph','alignment_graph.pkl'),
        os.path.join(OUTDIR,'log','alignment_graph_processing','weighted_alignment_graph.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'process_blastout.py'),
        walltime='86400',
        nodes='1',
        ppn='4'
    resources:
        mem = '4GB'
    threads:
        4
    log:
        o=os.path.join(OUTDIR,'qsub','alignment_graph_processing','weighted_alignment_graph.out'),
        e=os.path.join(OUTDIR,'qsub','alignment_graph_processing','weighted_alignment_graph.err'),
    conda:
        'vamb_n2v_master'
        #'envs/vamb.yaml'
    shell:
        """
        python {params.path} --blastout {input[0]} --out {output[0]}
        touch {output[1]}
        """


## 4. Merge assembly graphs an alignment graph into a unified graph
rule create_assembly_alignment_graph:
    input:
        alignment_graph_file = os.path.join(OUTDIR,'tmp','alignment_graph','alignment_graph.pkl'),
        assembly_graph_files =expand(os.path.join(OUTDIR,'tmp','assembly_graphs','{sample}.pkl'),sample=[sample for sample in SAMPLES.split()]),
        weighted_alignment_graph_finished_log = os.path.join(OUTDIR,'log','alignment_graph_processing','weighted_alignment_graph.finished'),
        weighted_assembly_graphs_all_samples_finished_log = os.path.join(OUTDIR,'log','assembly_graph_processing','weighted_assembly_graphs_all_samples.finished')
    output:
        os.path.join(OUTDIR,'tmp','assembly_alignment_graph.pkl'),
        os.path.join(OUTDIR,'log','create_assembly_alignment_graph.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'merge_assembly_alignment_graphs.py'),
        walltime='86400',
        nodes='1',
        ppn='4'
    resources:
        mem = '4GB'
    threads:
        4
    log:
        o=os.path.join(OUTDIR,'qsub','create_assembly_alignment_graph.out'),
        e=os.path.join(OUTDIR,'qsub','create_assembly_alignment_graph.err')
    conda:
        'vamb_n2v_master'
        #'envs/vamb.yaml'
    shell:
        """
        python {params.path} --graph_alignment {input.alignment_graph_file}  --graphs_assembly {input.assembly_graph_files} --out {output[0]} 
        touch {output[1]}
        """
        



## 3. Run n2v on the per sample assembly graphs
rule n2v_assembly_alignment_graph:
    input:
        os.path.join(OUTDIR,'tmp','assembly_alignment_graph.pkl'),
        os.path.join(OUTDIR,'log','create_assembly_alignment_graph.finished')
    output:
        directory(os.path.join(OUTDIR,'tmp','n2v','assembly_alignment_graph_embeddings')),
        os.path.join(OUTDIR,'tmp','n2v','assembly_alignment_graph_embeddings','embeddings.npz'),
        os.path.join(OUTDIR,'tmp','n2v','assembly_alignment_graph_embeddings','contigs_embedded.txt'),
        os.path.join(OUTDIR,'log','n2v','n2v_assembly_alignment_graph.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'fastnode2vec_args.py'),
        walltime='86400',
        nodes='1',
        ppn='8'
    resources:
        mem = '20GB'
    threads:
        8
    log:
        o=os.path.join(OUTDIR,'qsub','n2v','n2v_assembly_alignment_graph.out'),
        e=os.path.join(OUTDIR,'qsub','n2v','n2v_assembly_alignment_graph.err')
    conda:
        'vamb_n2v_master'
        #'envs/vamb.yaml'
    shell:
        """
        python {params.path} -G {input[0]} --ed {N2V_ED} --nw {N2V_NW} --ws {N2V_WS} --wl {N2V_WL}\
         -p {N2V_P} -q {N2V_Q} --outdirembs {output[0]} --normE {N2V_NZ} --contignames {CONTIGNAMES_FILE}
        touch {output[3]}
        """



## 4. Extract hoods from assembly graphs n2v embeddings per sample 
rule extract_neighs_from_n2v_embeddings:
    input:
        os.path.join(OUTDIR,'tmp','n2v','assembly_alignment_graph_embeddings','embeddings.npz'),
        os.path.join(OUTDIR,'tmp','n2v','assembly_alignment_graph_embeddings','contigs_embedded.txt'),
        os.path.join(OUTDIR,'log','n2v','n2v_assembly_alignment_graph.finished'),
        os.path.join(OUTDIR,'tmp','assembly_alignment_graph.pkl')

    output:
        directory(os.path.join(OUTDIR,'tmp','neighs')),
        os.path.join(OUTDIR,'tmp','neighs','neighs_object_r_%s.npz'%NEIGHS_R),
        os.path.join(OUTDIR,'log','neighs','extract_neighs_from_n2v_embeddings.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'embeddings_to_neighs.py'),
        walltime='86400',
        nodes='1',
        ppn='8'
    resources:
        mem = '50GB'
    threads:
        8
    log:
        o=os.path.join(OUTDIR,'qsub','neighs','extract_neighs_from_n2v_embeddings.out'),
        e=os.path.join(OUTDIR,'qsub','neighs','extract_neighs_from_n2v_embeddings.err')
    conda:
        'vamb_n2v_master'
    shell:
        """
        python {params.path} --embs {input[0]} --contigs_embs {input[1]}\
         --contignames {CONTIGNAMES_FILE} -g {input[3]} -r {NEIGHS_R} --neighs_outdir {output[0]}
        touch {output[2]}
        """

## 5. Run vamb to merge the hoods
rule run_vamb_asymmetric:
    input:
        os.path.join(OUTDIR,'log','neighs','extract_neighs_from_n2v_embeddings.finished'),
        CONTIGS_FILE,
        ABUNDANCE_FILE,
        os.path.join(OUTDIR,'tmp','neighs','neighs_object_r_%s.npz'%NEIGHS_R)#,
        #os.path.join(OUTDIR,'tmp','neighs','neighs_mask_r_%s.npz'%NEIGHS_R)

    output:
        #directory(os.path.join(OUTDIR,'vamb_asymmetric')),
        #os.path.join(OUTDIR,'vamb_asymmetric','vae_clusters_within_radius_complete_unsplit.tsv'),
        #os.path.join(OUTDIR,'vamb_asymmetric','vae_clusters_unsplit.tsv'),
        os.path.join(OUTDIR,'log/run_vamb_asymmetric.finished')
    params:
        walltime='86400',
        nodes='1',
        ppn=PLAMB_PPN,
        cuda='--cuda' if CUDA else ''
    resources:
        mem=PLAMB_MEM
    threads:
        int(plamb_threads)
    conda:
        'vamb_n2v_master' 
    log:
        o=os.path.join(OUTDIR,'qsub','run_vamb_asymmetric.out'),
        e=os.path.join(OUTDIR,'qsub','run_vamb_asymmetric.err')
    shell:
        """
        module unload gcc/13.2.0
        module unload gcc/12.2.0
        module load gcc/13.2.0
        {PLAMB_PRELOAD}
        vamb bin vae_asy --outdir {OUTDIR}/vamb_asymmetric --fasta {input[1]} -p {threads} --rpkm {input[2]}\
        --seed 1 --neighs {input[3]}  -m {MIN_CONTIG_LEN} {PLAMB_PARAMS}\
         {params.cuda}  
        touch {output}
        """


rule run_geNomad:
    input:
        CONTIGS_FILE
    output:
        directory(os.path.join(OUTDIR,'tmp','geNomad')),
        os.path.join(OUTDIR,'log/run_geNomad.finished'),
        os.path.join(OUTDIR,'tmp','geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv')
    params:
        db_geNomad="/home/projects/ku_00197/data/ptracker/dbs/genomad_db",
        walltime='86400',
        nodes='1',
        ppn='16'
    resources:
        mem='50GB'
    threads:
        16
    log:
        o = os.path.join(OUTDIR,'qsub/run_geNomad.o'),
        e = os.path.join(OUTDIR,'qsub/run_geNomad.e')
    conda:
        #'envs/blast.yaml'
        'vamb_n2v_master'
    shell:
        """
        module load mmseqs2/release_15-6f452
        module load aragorn/1.2.36
        genomad end-to-end --cleanup {input} {output[0]}   {params.db_geNomad} --threads {threads}
        touch {output[1]}
        """


## 7. Classify bins/clusters into plasmid/organism/virus bins/clusters
rule classify_bins_with_geNomad:
    input:
        os.path.join(OUTDIR,'tmp','geNomad','contigs.flt_aggregated_classification','contigs.flt_aggregated_classification.tsv'),
        os.path.join(OUTDIR,'log/run_vamb_asymmetric.finished'),
        os.path.join(OUTDIR,'log/run_geNomad.finished')
    output:
        os.path.join(OUTDIR,'vamb_asymmetric','vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv'),
        os.path.join(OUTDIR,'log','classify_bins_with_geNomad.finished')
    params:
        path = os.path.join(SNAKEDIR, 'src', 'classify_bins_with_geNomad.py'),
        walltime='86400',
        nodes='1',
        ppn=1,
    resources:
        mem="1GB"
    threads:
        1
    conda:
        'vamb_n2v_master' 
    log:
        o=os.path.join(OUTDIR,'qsub','classify_bins_with_geNomad.out'),
        e=os.path.join(OUTDIR,'qsub','classify_bins_with_geNomad.err')
    shell:
        """
        python {params.path} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit.tsv\
         --dflt_cls {OUTDIR}/vamb_asymmetric/vae_clusters_unsplit.tsv --scores {input[0]} --outp {output[0]} --lengths {OUTDIR}/vamb_asymmetric/lengths.npz --contignames {OUTDIR}/vamb_asymmetric/contignames
        touch {output[1]}
        """
    


rule rename_clusters_and_hoods:
    input:
        os.path.join(OUTDIR,'log','classify_bins_with_geNomad.finished')
    output:
        os.path.join(OUTDIR,'log/rename_clusters_and_hoods.finished')
    params:
        path_sample = os.path.join(SNAKEDIR, 'src', 'replace_samples.py'),
        path_cluster = os.path.join(SNAKEDIR, 'src', 'rename_clusters.py'),
        walltime='86400',
        nodes='1',
        ppn='1'

    resources:
        mem="1GB"
    threads:
        1
    conda:
        'vamb_n2v_master' 
    log:
        o=os.path.join(OUTDIR,'qsub','rename_clusters_and_hoods.out'),
        e=os.path.join(OUTDIR,'qsub','rename_clusters_and_hoods.err')
    shell:
        """
        python {params.path_cluster} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit_renamed.tsv --header
        python {params.path_cluster} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_unsplit.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_unsplit_renamed.tsv --header
        python {params.path_cluster} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_unsplit.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/vamb_asymmetric/vae_clusters_unsplit_renamed.tsv  --header 
        
        python {params.path_cluster} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/vamb_asymmetric/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids_renamed.tsv --header
        python {params.path_cluster} --clusters {OUTDIR}/vamb_asymmetric/vae_clusters_unsplit_geNomadplasclustercontigs_extracted.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/vamb_asymmetric/vae_clusters_unsplit_geNomadplasclustercontigs_extracted_renamed.tsv  --header 
        
        python {params.path_cluster} --clusters {OUTDIR}/tmp/neighs/hoods_clusters_r_0.05.tsv --renaming_file {RENAMING_FILE} --out {OUTDIR}/tmp/neighs/hoods_clusters_r_0.05_renamed.tsv --header
        touch {output}
        """

