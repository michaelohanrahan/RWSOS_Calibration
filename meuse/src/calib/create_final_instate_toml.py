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
                  ST_key,
                  soilthickness,
                  start,
                  end,
                  log):
    
    # this is a toml file
    config = model.config
    # update the log 
    config["log"] = f"../0-log/instate_final/instates.txt"
    
    # update the reinit
    config["model"]["reinit"] = True
    
    # update the input path_static
    config["input"]["path_static"] = f"staticmaps.nc"
    
    # update the input path_forcing
    config["input"]["path_forcing"] = "forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc"
    
    # print(config.keys())
    # update the state path_output
    config["state"]["path_output"] = f'instate.nc'
    config["state"]["path_input"] = None

    # Update the starttime
    config["starttime"] = start
    
    # Update the endtime
    config["endtime"] = end
    
    #we dont need the output from this period, csv outputs of generic name compete 
    #stopping them from running in parallel
    config.pop("csv", None)
    # print(config.keys())
    
    # Write the config to a TOML file
    output_path = Path(root) / "instates" / f"wflow_sbm_getinstate_final.toml"
    with open(output_path, 'w', encoding='utf-8') as toml_file:
        toml.dump(config, toml_file)

    return model



if __name__ == "__main__":
    
    l = setup_logging('data/0-log', f'01-initial_instate_tomls_L{snakemake.params.level}.log')
    try:
        config_fn = snakemake.input.config_fn
        root= snakemake.params.root
        start = snakemake.params.starttime
        end = snakemake.params.endtime
        os.makedirs("instates", exist_ok=True)
        l.info(f"Params:\n config_fn: {config_fn}\n root: {root}\n ST_values: {ST_values}\n ST_key: {ST_key}\n level: {level}\n start: {start}\n end: {end}")
        l.info(f"Creating instate files for soil thickness values: {ST_values}")
    
    except NameError:
        os.chdir(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse')
        l.info("Running script in test mode")
        config_fn= Path("wflow_sbm.toml")   
        root= Path("data/3-input")
        ST_values =[999, 9999]
        ST_key = "wrongsoilthicknesskey"
        level = 0
        start = "2005-01-01T00:00:00"
        end = "2005-03-01T00:00:00"
        os.makedirs("instates", exist_ok=True)
    
    except Exception as e:
        l.error(f"An error occurred importing params: {e}")
        l.error(traceback.format_exc())
        raise e

    try:
        mod = read_model(config_fn, l)
        l.info(f'original mod.config.input: {mod.config["input"]}')
        change_config(model=mod,
                    root=root, 
                    start=start, 
                    end=end,
                    log=l)
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e
    
        
    
        
