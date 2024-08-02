#!/usr/bin/env julia

using BinBencherBackend

reference = ARGS[1]#"/home/projects/ku_00197/data/ptracker/tmp/outdir_refs_ef/ref_spades_$location.json"
only_organism = ARGS[2]#"false"
bins = ARGS[3]#"/home/people/lasdan/lasdan/default_vamb_lasse/Airways//vamb_n2v_master_default_261123/261123_vamb_n2v_master_default_1/vae_clusters_split.tsv"
assembly = ARGS[4]#"false"
level = 1

ref = Reference(reference)

if only_organism == "true"
  subset!(ref;
            genomes = is_organism
              )
end


bins = open(i -> Binning(i, ref), bins)

if assembly == "true"
  nc_strains = n_recovered(bins, 0.90, 0.95; level=level, assembly=true)
elseif assembly == "false"
  nc_strains = n_recovered(bins, 0.90, 0.95; level=level, assembly=false)
end
print(ARGS[3], "\t",only_organism, "\t",assembly,"\t", nc_strains, "\n")