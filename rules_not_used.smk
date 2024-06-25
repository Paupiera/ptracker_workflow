

# rulename = "download"
# rule download:
#     output:
#         #"data/sample_{key}/fastq/{id}_1.fastq", # Quick fix for now.. improve it later
#         #"data/sample_{key}/fastq/{id}_2.fastq",
#         fw = "data/reads/{id}_1" + config["fastq_end"], # Quick fix for now.. improve it later .. actually not worth it.. i guess
#         rv = "data/reads/{id}_2" + config["fastq_end"],
#         prefetched_info = temp("{id}"),
#     params:
#         id_download = "data/sample_{key}/fastq",
#     group: "spades"
#     benchmark: config["benchmark"]+rulename
#     log: e =  config["log.e"]+rulename, o =  config["log.o"]+rulename,
#     shell:
#         "bin/sratoolkit.3.0.10-ubuntu64/bin/prefetch {wildcards.id} " 
#         +LOG_CMD+ 
#         "bin/sratoolkit.3.0.10-ubuntu64/bin/fastq-dump {wildcards.id} -O {params.id_download} --split-3 " 
#         +LOG_CMD+
        #"rm -rf {wildcards.id}" # Remove the pre-fetched information after collecting the fastq files


##  START PAU SETUP
def get_config(name, default, regex):
    #res = config.get(name, default).strip()
    CONFIGFILE = "PauPipeline/config_Gastrointestinal_with_renaming_geNomad.json" # TODO delete and add info to the other configfile
    res = CONFIGFILE
    m = re.match(regex, res)
    if m is None:
        raise ValueError(
            f'Config option \'{name}\' is \'{res}\', but must conform to regex \'{regex}\'')
    return res

# rule split_reads:
#  input:
#        paired = "data/reads/{id}.fq.gz",
#  output:
#        fw = read_fw, # "data/reads/{id}_1" + config["fastq_end"],
#        rv = read_rv, #"data/reads/{id}_2" + config["fastq_end"],
#  threads: 8
#  shell:
#        """
#         zcat {input.paired} | \
#         paste - - - - - - - -  | tee >(cut -f 1-4 | tr "\t" "\n" | \
#         pigz --best --processes {threads} > {output.fw}) | \
#         cut -f 5-8 | tr "\t" "\n" | pigz --best --processes {threads} > {output.rv}
#         # https://gist.github.com/nathanhaigh/3521724
#        """

# rulename="genomad"
# rule genomad:
#     input: 
#         fasta = "data/sample_{key}/contigs.flt.fna.gz", 
#     output:
#         direc = "data/sample_{key}/genomad",
#         file = "data/sample_{key}/genomad/contigs.flt_aggregated_classification",
#     params:
#         db = "genomad_db", 
#     threads: config["genomad"]["threads"]
#     resources: walltime = config["genomad"]["walltime"], mem_gb = config["genomad"]["mem_gb"]
#     #conda:
#     #    "genomad"
#     benchmark: config["benchmark.key"]+rulename
#     log: e =  config["log.e.key"]+rulename, o =  config["log.o.key"]+rulename,
#     shell:
#         "genomad end-to-end {input.fasta} {output.direc} {params.db} "
#         "-t {threads}"
#         + LOG_CMD