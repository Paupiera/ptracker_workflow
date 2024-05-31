
shell.prefix("")
rule test:
    output: "delme.txt"
    shell: 
        """
        touch {output}
        """


# import pandas as pd
# configfile: "config/config.yaml"
# from collections import UserDict




# def populate_dict_of_lists(dict_to_pop, key, id):
#     if key not in dict_to_pop:
#         dict_to_pop[key] = [id]
#     else:
#         dict_to_pop[key].append(id)

# df = pd.read_csv("./config/accessions.txt", sep="\s+", comment="#")
# sample_id = {}
# sample_path = {}
# for sample, id in zip(df.SAMPLE, df.ID):
#     id = str(id)
#     sample = str(sample)

#     populate_dict_of_lists(sample_id, sample, id)

# print(sample_id)
# print(sample_path)






# def get_config(name, default, regex):
#     #res = config.get(name, default).strip()
#     CONFIGFILE = "PauPipeline/config_Gastrointestinal_with_renaming_geNomad.json" # TODO delete and add info to the other configfile
#     res = CONFIGFILE
#     m = re.match(regex, res)
#     if m is None:
#         raise ValueError(
#             f'Config option \'{name}\' is \'{res}\', but must conform to regex \'{regex}\'')
#     return res



