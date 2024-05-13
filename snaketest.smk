configfile: "config/config.yaml"
print(config)

class config_dict(dict):
    pass
    




def get_config(name, default, regex):
    #res = config.get(name, default).strip()
    CONFIGFILE = "PauPipeline/config_Gastrointestinal_with_renaming_geNomad.json" # TODO delete and add info to the other configfile
    res = CONFIGFILE
    m = re.match(regex, res)
    if m is None:
        raise ValueError(
            f'Config option \'{name}\' is \'{res}\', but must conform to regex \'{regex}\'')
    return res



