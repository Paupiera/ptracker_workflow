
import sys
import pandas as pd
from collections import UserDict



def expand_dir(torep, dic):
    allcomb = []
    for key, val in dic.items():
        tmptorep = torep
        tmptorep = tmptorep.replace("[key]", key)
        for value in val:
            allcomb.append(tmptorep.replace("[value]", value))
    return allcomb

def populate_dict_of_lists(dict_to_pop, key, id):
    if key not in dict_to_pop:
        dict_to_pop[key] = [id]
    else:
        dict_to_pop[key].append(id)

class config_dict(UserDict):
    """
    Dict which also support a method to get the value, or return a default value if it's none
    """
    def __init__(self, init_dict: dict):
        """
        Convert a dict to a config_dict, recursively
        """

        self_data = dict()
        for key, value in init_dict.items():
            if type(value) == dict:
                self_data[key] = config_dict(value)
            else:
                self_data[key] = value
        self.data = self_data


    def get_item_or_default(self, key, default_value):
        if self.data[key] == None:
            return default_value
        return self.data[key]

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
