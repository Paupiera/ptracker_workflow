conda envs er funky -- suprise
try with module snakemake -- in the future all conda should be called through snakemake, 
and snakemake can then be it's own env
- does not work as the module snakemake messes everything up
- Try with conda env through snakemake and pray to god
 

# TODO megahit istedet for spades?
# TODO Downloading downloas 3 files, _1 and _2 and a second - what is it.. and should it be used (likey non-paired reads something)
# TODO Compress non-interleved fastq files


df = pd.read_csv("./config/test_accessions.txt", sep="\s+", comment="#")
sample_id = {}
sample_path = {}
for sample, id, path in zip(df.SAMPLE, df.ID, df.PATH):
    id = str(id)
    path = str(path)
    sample = str(sample)

    populate_dict_of_lists(sample_id, sample, id)
    populate_dict_of_lists(sample_path, sample, path)