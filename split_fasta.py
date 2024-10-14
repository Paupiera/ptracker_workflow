import click

def readfasta(filename):
    '''Reading in several fasta files
       yield each pair of header and dna seperate
    '''
    # Try to open the file, error if file does not exsist
    try: file = open(filename, "r")
    except FileNotFoundError as errormessage:
        sys.exit(f"The file '{filename}' could not be found, error: {errormessage}")

    # Extract the dna, and the headers.
    oldheader = "FIRST HEADER"
    for line in file:
        line = line.strip()
        #if line is header yield dna and header except FIRST HEADER
        if line.startswith(">"):    
            newheader = line
            if oldheader != "FIRST HEADER":
                yield dna, oldheader
            dna = ""
            oldheader = newheader
        else:
            dna += line
    # Yield the last header and dna
    yield dna, oldheader
    file.close()

import click
@click.command()
@click.option('--fasta_all_contigs', type=click.Path(exists=True), required = True)
@click.option('--clusterfile_plasmid', type=click.Path(exists=True), required = True)
@click.option('--clusterfile_chromosome', type=click.Path(exists=True), required = True)
def split_fasta(fasta_all_contigs, clusterfile_plasmid, clusterfile_chromosome):
    import pandas as pd                                                              
    import sys
    import os

    # Example files on ESRUM
    # clusterfile_plasmid = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv" 
    # clusterfile_chromosome = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv" 
    # fasta_all_contigs = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/contigs_darwin_sub.fna"                           

    df_chromosome = pd.read_table(clusterfile_chromosome, sep="\t")                                             
    df_plasmid = pd.read_table(clusterfile_plasmid, sep="\t")                                             

    df_plasmid = df_plasmid.assign(contigname = ">" + df_plasmid.contigname)
    df_chromosome = df_chromosome.assign(contigname = ">" + df_chromosome.contigname)

# Get the sum of the lengths of each contig in each bin
    df_chromosome = df_chromosome.assign(contig_length = df_chromosome.contigname.str.extract(r"length_(\d+)").astype("int"))
    df_chromosome["cluster_length"] = df_chromosome.groupby("clustername")["contig_length"].transform("sum")

# Create 3 sets to filter fasta files later, containing different contigs
    t = 200_000
    chromosomes_above_t = set(df_chromosome.contigname[df_chromosome.cluster_length >= t])
    chromosomes_below_t = set(df_chromosome.contigname[df_chromosome.cluster_length < t])
    plasmids = set(df_plasmid.contigname)

    plasmid_dir = "candidate_plasmids"
    os.mkdir("candidate_plasmids")
    chromosome_dir = "candidate_chromosomes"
    os.mkdir("candidate_chromosomes")

# Write to the split fastafiles
    with open(f'{plasmid_dir}/plasmids.fna', 'w') as plasmid_file, \
        open(f'{chromosome_dir}/chromosomes_candidate.fna', 'w') as chromosome_above_t, \
        open('chromosome_contigs_filtered_away.fna', 'w') as chromosome_below_t:
        for dna, header in readfasta(fasta_all_contigs ):
            if header in plasmids:
                    print(header, file=plasmid_file)
                    print(dna, file=plasmid_file)
            elif header in chromosomes_above_t:
                    print(header, file=chromosomes_above_t)
                    print(dna, file=chromosomes_above_t)
            elif header in chromosomes_below_t:
                    print(header, file=chromosomes_below_t)
                    print(dna, file=chromosomes_below_t)



if __name__ == "__main__":
    split_fasta()
