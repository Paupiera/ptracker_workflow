
import sys
import pandas as pd



def expand_dir(torep, dic):
    allcomb = []
    for key, val in dic.items():
        tmptorep = torep
        tmptorep = tmptorep.replace("{key}", key)
        for value in val:
            

            allcomb.append(tmptorep.replace("{value}", value))
    return allcomb

# file = "./config/accessions.txt"
# with open(file, "r") as f:
#     out = [line.strip() for line in f.readlines() if line.strip() != ""]
#     print(out)

# df = pd.read_csv("./config/accessions.txt", sep="\s+")
# sample_id = {}
# for sample,id in zip(df.SAMPLE, df.ID):
#     if sample not in sample_id:
#         sample_id[sample] = [id]
#     else:
#         sample_id[sample].append(id)
#
# print(sample_id)
#
#
# df.to_csv(sys.stdout, sep="\t")
