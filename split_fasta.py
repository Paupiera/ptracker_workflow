# import click
from pathlib import Path

from numpy import single
import rich_click as click
import gzip

def readfasta(filename):
    '''Reading in several fasta files
       yield each pair of header and dna seperate

    '''
    if filename.endswith('.gz'):
        import gzip
        my_open = gzip.open
    else:
        my_open = open

    # Try to open the file, error if file does not exsist
    try: file = my_open(filename, "rt")
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

def print_fasta_by_group(df, togroupby, dir, count, header2dna):
    dir = Path(dir)
    output_path = dir / count.__repr__()

    with open(output_path, "w") as f:
        for contigname in df['contigname']:
            print(contigname, file=f)
            print(header2dna[contigname], file=f)

    count.add_one()


class Count():
    count = 0
    def add_one(self):
        self.count += 1
    def __repr__(self):
        return str(self.count)


def write_clusters(df, dir, header2dna):
    # Write each cluster to a dir 
    count = Count()
    togroupby = "clustername"
    (df
        .groupby(togroupby)
        .agg(print_fasta_by_group, togroupby, dir, count, header2dna)
        )

@click.command()
@click.option('--fasta_all_contigs', type=click.Path(exists=True), required = True, help="Fasta file containing all contigs to sort through")
@click.option('--clusterfile_plasmid', type=click.Path(exists=True), required = True, help="Clusterfile in vambs output format of all plasmids")
@click.option('--clusterfile_chromosome', type=click.Path(exists=True), required = True, help="Clusterfile in vambs output format of all non-plasmids")
@click.option('--checkm_output', type=click.Path(exists=False) ,help="Output the format in checkm2 format to the dir given as argument")
@click.option('--single_file_output', type=click.Path(exists=False) ,help="Output the format in singlefile output")
def split_fasta(fasta_all_contigs, clusterfile_plasmid, clusterfile_chromosome, checkm_output, single_file_output):
    import pandas as pd                                                              
    import sys
    import os

    # Example files on ESRUM
    # clusterfile_plasmid = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv" 
    # clusterfile_chromosome = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/vae_clusters_within_radius_with_looners_complete_unsplit_candidate_plasmids.tsv" 
    # fasta_all_contigs = "/home/bxc755/rasmussen/scratch/ptracker/plasmid_graph/data/contigs_darwin_sub.fna"                           

    df_plasmid = pd.read_table(clusterfile_plasmid, sep="\t")                                             
    # TODO change this to filter out the ones containing only plasmids
    plasmids = set(df_plasmid.contigname)
    df_chromosome = pd.read_table(clusterfile_chromosome, sep="\t")                                             

    df_plasmid = df_plasmid.assign(contigname = ">" + df_plasmid.contigname)
    df_chromosome = df_chromosome.assign(contigname = ">" + df_chromosome.contigname)

    # Get the sum of the lengths of each contig in each bin
    df_chromosome = df_chromosome.assign(contig_length = df_chromosome.contigname.str.extract(r"length_(\d+)").astype("int"))
    df_chromosome["cluster_length"] = df_chromosome.groupby("clustername")["contig_length"].transform("sum")

    if single_file_output != None:
        plasmid_dir = "candidate_plasmids"
        Path(plasmid_dir).mkdir(exist_ok=True)
        chromosome_dir = "candidate_chromosomes"
        Path(chromosome_dir).mkdir(exist_ok=True)

        # Create 3 sets to filter fasta files later, containing different contigs
        t = 200_000
        chromosomes_above_t = set(df_chromosome.contigname[df_chromosome.cluster_length >= t])
        chromosomes_below_t = set(df_chromosome.contigname[df_chromosome.cluster_length < t])

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

    
    if checkm_output != None:

        # Make various paths
        checkm_output = Path(checkm_output)
        checkm_output.mkdir(exist_ok=True)
        chromosome_dir = checkm_output / "candidate_chromosomes"
        non_candidate_chromosome_dir = checkm_output / "non_candidate_chromosomes"
        plasmid_dir = checkm_output / "candidate_plasmids"
        plasmid_dir.mkdir()
        chromosome_dir.mkdir()
        non_candidate_chromosome_dir.mkdir()

        # Read in the fasta files
        header2dna = {header: dna for dna, header in readfasta(fasta_all_contigs)}

        # Write each cluster to a dir for the plasmid bins
        write_clusters(df_plasmid, plasmid_dir, header2dna)

        # and for the chromsomes
        t = 200_000
        chromosomes_above_t = df_chromosome[df_chromosome.cluster_length >= t]
        print(chromosomes_above_t)
        chromosomes_below_t = df_chromosome[df_chromosome.cluster_length < t]
        write_clusters(chromosomes_above_t, chromosome_dir, header2dna)
        write_clusters(chromosomes_below_t, non_candidate_chromosome_dir, header2dna)


if __name__ == "__main__":
    split_fasta()
