# TODO path to files should not be thorugh symlinks -- access for other people
# TODO add genomad til mixed
# TODO megahit istedet for spades?
# TODO Downloading downloas 3 files, _1 and _2 and a second - what is it.. and should it be used (likey non-paired reads something)
configfile: "config/config.yaml"

from utils import expand_dir
import pandas as pd

# Constants
LOG_CMD = " > {log.o} 2> {log.e};"

# Read in the sample data
df = pd.read_csv("./config/accessions.txt", sep="\s+", comment="#")
sample_id = {}
for sample,id in zip(df.SAMPLE, df.ID):
    if sample not in sample_id:
        sample_id[sample] = [id]
    else:
        sample_id[sample].append(id)

# print( expand_dir("data/sample_{key}/mapped/{value}.bam.sort", sample_id))
rule all:
    input:
        expand_dir("data/sample_{key}/mapped/{value}.bam.sort", sample_id),
        expand("data/sample_{key}/genomad/contigs.flt_aggregated_classification", key=sample_id.keys()),
        # expand("data/sample_{key}/mapped/{value}.sort.bam", key=sample_id.keys())

#rulename = "download"
#rule download:
#    output:
#        "data/sample_{key}/fastq/{id}_1.fastq", # Quick fix for now.. improve it later
#        "data/sample_{key}/fastq/{id}_2.fastq",
#    params:
#        id_download = "data/sample_{key}/fastq",
#    group: "spades"
#    benchmark: config["benchmark"]+rulename
#    log: e =  config["log.e"]+rulename, o =  config["log.o"]+rulename,
#    shell:
#        "/home/people/lasdan/ptracker/sratoolkit.3.0.10-ubuntu64/bin/prefetch {wildcards.id} " 
#        +LOG_CMD+ 
#        "/home/people/lasdan/ptracker/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump {wildcards.id} -O {params.id_download} --split-3 " 
#        +LOG_CMD+
#        "rm -rf {wildcards.id}" # Remove the pre-fetched information after collecting the fastq files

rulename = "spades"
rule spades:
    input:
        #fw = "data/sample_{key}/fastq/{id}_1.fastq",    |
        #rv = "data/sample_{key}/fastq/{id}_2.fastq",    | -- For when rule download is uncommented
        fw = "data/reads/{id}_1.fastq", # Quick fix for now.. improve it later .. actually not worth it.. i guess
        rv = "data/reads/{id}_2.fastq",
    output:
        outdir = directory("data/sample_{key}/spades_{id}"),
        outfile = "data/sample_{key}/spades_{id}/contigs.fasta",
    group: "spades"
    threads: config["spades"]["threads"]
    resources: walltime = config["spades"]["walltime"], mem_gb = config["spades"]["mem_gb"]
    benchmark: config["benchmark"]+rulename
    log: e =  config["log.e"]+rulename, o =  config["log.o"]+rulename,
    shell:
        #"rm -rf {output.outdir};"
        "bin/SPAdes-3.15.4-Linux/bin/metaplasmidspades.py "
        "-o {output.outdir} -1 {input.fw} -2 {input.rv} " 
        "-t {threads} " # No reason to set mem, as spades just quits if it uses above the threshold
        +LOG_CMD

rulename="cat_contigs"
rule cat_contigs:
    input:
        expand_dir("data/sample_{key}/spades_{value}/contigs.fasta", sample_id)
    output:
        "data/sample_{key}/contigs.flt.fna.gz"
    #conda: 
    #    "vambworks"
    group: "spades"
    benchmark: config["benchmark.key"]+rulename
    log: e =  config["log.e.key"]+rulename, o =  config["log.o.key"]+rulename,
    shell: 
        "python bin/vamb/src/concatenate.py {output} {input} -m 2000" 
        +LOG_CMD

# Index resulting contig-file with minimap2
rulename="index"
rule index:
    input:
        contigs = "data/sample_{key}/contigs.flt.fna.gz"
    output:
        mmi = "data/sample_{key}/contigs.flt.mmi"
    group: "spades"
    #envmodules:
    #    'tools',
    #    'minimap2/2.6'
    benchmark: config["benchmark.key"]+rulename
    log: e =  config["log.e.key"]+rulename, o =  config["log.o.key"]+rulename,
    shell:
        "minimap2 -d {output} {input} "
        +LOG_CMD


# This rule creates a SAM header from a FASTA file.
# We need it because minimap2 for truly unknowable reasons will write
# SAM headers INTERSPERSED in the output SAM file, making it unparseable.
# To work around this mind-boggling bug, we remove all header lines from
# minimap2's SAM output by grepping, then re-add the header created in this
# rule.
rulename="dict"
rule dict:
    input:
        contigs = "data/sample_{key}/contigs.flt.fna.gz",
    output:
        "data/sample_{key}/contigs.flt.dict"
    group: "spades"
    #envmodules:
    #    'tools',
    #    'samtools/1.17',
    benchmark: config["benchmark.key"]+rulename
    # log: e =  config["log.e.key"]+rulename, o =  config["log.o.key"]+rulename,
    shell:
        "samtools dict {input.contigs} | cut -f1-3 > {output}"

# Generate bam files 
LONG_READS = False
rulename="minimap"
rule minimap:
    input:
        #fw = "data/sample_{key}/fastq/{id}_1.fastq", |
        #rv = "data/sample_{key}/fastq/{id}_2.fastq", | - For when downloading fastq files
        fw = "data/reads/{id}_1.fastq", # Quick fix for now.. improve it later .. actually not worth it.. i guess
        rv = "data/reads/{id}_2.fastq",

        mmi ="data/sample_{key}/contigs.flt.mmi",
        dict = "data/sample_{key}/contigs.flt.dict"
    output:
        bam = temp("data/sample_{key}/mapped/{id}.bam")
    group: "minimap"
    params:
        long_or_short_read = 'map-pb -L' if LONG_READS else 'sr',
    #envmodules:
    #    'tools',
    #    'samtools/1.17',
    #    'minimap2/2.6'
    threads: config["minimap"]["threads"]
    resources: walltime = config["minimap"]["walltime"], mem_gb = config["minimap"]["mem_gb"]
    benchmark: config["benchmark"]+rulename
    # log: e =  config["log.e"]+rulename, o =  config["log.o"]+rulename,
    shell:
        # See comment over rule "dict" to understand what happens here
        "minimap2 -ax {params.long_or_short_read} {input.mmi} {input.fw} {input.rv} "
        " | grep -v '^@'"
        " | cat {input.dict} - "
        " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
        " > {output.bam}"

 # Sort bam files
rulename="sort"
rule sort:
    input:
        "data/sample_{key}/mapped/{id}.bam",
    output:
        "data/sample_{key}/mapped/{id}.bam.sort",
    group: "minimap"
    #envmodules:
    #    'tools',
    #    'samtools/1.17',
    benchmark: config["benchmark"]+rulename
    # log: e =  config["log.e"]+rulename, o =  config["log.o"]+rulename,
    shell:
        "samtools sort {input} -o {output} "

rulename="genomad"
rule genomad:
    input: 
        fasta = "data/sample_{key}/contigs.flt.fna.gz", 
    output:
        direc = "data/sample_{key}/genomad",
        file = "data/sample_{key}/genomad/contigs.flt_aggregated_classification",
    params:
        db = "genomad_db", 
    threads: config["genomad"]["threads"]
    resources: walltime = config["genomad"]["walltime"], mem_gb = config["genomad"]["mem_gb"]
    #conda:
    #    "genomad"
    benchmark: config["benchmark.key"]+rulename
    log: e =  config["log.e.key"]+rulename, o =  config["log.o.key"]+rulename,
    shell:
        "genomad end-to-end {input.fasta} {output.direc} {params.db} "
        "-t {threads}"
        + LOG_CMD






