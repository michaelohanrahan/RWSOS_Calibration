# import snakemake
from pathlib import Path
from hydromt_wflow import WflowModel
from hydromt.log import setuplog
import traceback
from setuplog import setup_logging
import os
import toml

def read_model(config_fn, log):
    # read the model configuration
    model = WflowModel(root=root, mode="r", config_fn=config_fn,  logger=log)
    model.read_config()
    # model.read_grid()
    return model

def change_config(model, 
                  root,
                  level, 
                  start,
                  end,
                  log,
                  outcfg):
    # this is a toml file
    config = model.config
    l.info(F"setting up instate for level {level}")
    # update the log 
    config["log"] = f"../../0-log/instate_L{level}/instates_level{level}.txt"
    
    # casename
    config["casename"] = f"instates_level{level}"
    
    # update the reinit
    config["model"]["reinit"] = True
    
    # update the input path_static
    config["input"]["path_static"] = f"../staticmaps.nc"
    
    # update the input path_forcing
    config["input"]["path_forcing"] = "../forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc"
    
    # update the state path_output
    config["state"]["path_output"] = f'instate_level_{level}.nc'
    config["state"]["path_input"] = None

    # Update the starttime
    config["starttime"] = start
    
    # Update the endtime
    config["endtime"] = end
    
    #we dont need the output from this period, csv outputs of generic name compete 
    #stopping them from running in parallel
    if "csv" in config:
        config.pop("csv", None)
    
    # Write the config to a TOML file
    output_path = outcfg
    with open(output_path, 'w', encoding='utf-8') as toml_file:
        toml.dump(config, toml_file)
    l.info(f"Updated config saved to {output_path}")
    return model



if __name__ == "__main__":
    
    l = setup_logging('data/0-log', f'01-initial_instate_tomls_L{snakemake.params.level}.log')
    try:
        config_fn = snakemake.input.config_fn
        root= snakemake.params.root
        level = snakemake.params.level
        start = snakemake.params.starttime
        end = snakemake.params.endtime
        os.makedirs("instates", exist_ok=True)
        
    
    except Exception as e:
        l.error(f"An error occurred importing params: {e}")
        l.error(traceback.format_exc())
        raise e

    try:
        mod = read_model(config_fn, l)
        change_config(model=mod,
                    root=root, 
                    level=level, 
                    start=start, 
                    end=end,
                    log=l,
                    outcfg=snakemake.output.cfg)
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e
    
        
    
        
