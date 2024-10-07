# import snakemake
from pathlib import Path
from hydromt_wflow import WflowModel
from hydromt.log import setuplog
import traceback
from setuplog import setup_logging
import os
import toml

def read_model(root, config_fn):
    # read the model configuration
    model = WflowModel(root=root, mode="r", config_fn=config_fn)
    model.read_config()
    # model.read_grid()
    return model

def change_config(model,
                  start,
                  end,
                  forcing_path,
                  topx,
                  l,
                  outcfg):
    # this is a toml file
    config = model.config
    l.info(F"setting up instate for {topx}")
    # update the log 
    config["log"] = f"../../0-log/instate_{topx}/instates_{topx}.txt"
    
    # casename
    config["casename"] = f"instates_{topx}"
    
    # update the reinit
    config["model"]["reinit"] = True
    
    # update the input path_static
    config["input"]["path_static"] = f"../staticmaps.nc"
    
    # update the input path_forcing
    config["input"]["path_forcing"] = forcing_path.as_posix()
    
    # update the state path_output
    config["state"]["path_output"] = f'instate_final.nc'
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
    
    l = setup_logging('data/0-log', f'01-initial_instate_tomls_{snakemake.params.topx}.log')
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        wflow_mod = read_model(mod.params.root, mod.input.config_fn)
        change_config(model=wflow_mod, 
                    start=mod.params.starttime, 
                    end=mod.params.endtime,
                    topx=mod.params.topx,
                    forcing_path=mod.params.forcing_path,
                    l=l,
                    outcfg=mod.output.cfg)
        
        
        
    
        
