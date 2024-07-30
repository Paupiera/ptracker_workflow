import re
import sys
file = sys.argv[1]

rulename = []
cluster_id = []
with open(file, "r") as f:
    for line in f:
        match_name = re.findall(r"Error executing rule (\w+) on cluster", line)
        match_cluster_id = re.findall(r"Submitted batch job (\d+)", line)

        if len(match_name) > 0:
            rulename.append(match_name[0])
        if len(match_cluster_id) > 0:
            cluster_id.append(match_cluster_id[0])
        
for name, id in zip(rulename, cluster_id):
    log_file = f"snakemake_output/{name}.{id}"
    print(log_file)



