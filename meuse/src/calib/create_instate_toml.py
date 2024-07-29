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

# def change_config(model, 
#                   root,
#                   level, 
#                   ST_key,
#                   soilthickness,
#                   start,
#                   end,
#                   log):
    
#     # Working
#     model.read_config()
#     model.set_config("log", f"../../0-log/instates_level{level}_ST{str(soilthickness).replace('.', '')}.log")
#     model.set_config("reinit", True)
#     model.set_config('input', {'path_static':"../staticmaps/staticmaps.nc"})
#     l.info(f"static file: {model.config['input']['path_static']}")
    
#     model.set_config('input', 'path_forcing', "../forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc")
#     l.info(f"forcing file: {model.config['input']['path_forcing']}")
    
#     model.set_config("state", {"path_output":f'instate_level{level}_ST{str(soilthickness).replace(".", "")}.nc'})
#     l.info(f"state file: {model.config['state']['path_output']}")
    
#     #notworking
    
#     #declare it as a nc variable in the 
#     model.set_config('input', 'vertical', 'soilthickness', 'netcdf', 'variable', 'name', ST_key)
#     l.info(f'nc variable name: {model.config["input"]["vertical"]["soilthickness"]["netcdf"]["variable"]["name"]}')
    
#     model.set_config('input', 'vertical', 'soilminthickness', 'netcdf', 'variable', {'name': ST_key}) 
#     l.info(f'nc variable name: {model.config["input"]["vertical"]["soilminthickness"]["netcdf"]["variable"]["name"]}')
       
#     model.set_config('input','vertical','soilthickness',{"scale":soilthickness})
#     l.info(f"soilthickness: {model.config['input']['vertical']['soilthickness']}")
    
#     model.set_config('input','vertical','soilminthickness',{"scale":soilthickness})
#     l.info(f"soilminthickness: {model.config['input']['vertical']['soilminthickness']}")
    
#     model.set_config("starttime", start)
#     l.info(f"start time: {model.config['starttime']}")
    
#     model.set_config("endtime", end)     
#     l.info(f"end time: {model.config['endtime']}")
    
    
    
#     model.set_root(Path(root, "instates"), mode="w+")
    
#     l.info(f"Writing instate file for soil thickness: {soilthickness}")
#     model.write_config(f"wflow_sbm_getinstate_level{level}_ST{str(soilthickness).replace('.', '')}.toml")
#     return model

def change_config(model, 
                  root,
                  level, 
                  ST_key,
                  soilthickness,
                  start,
                  end,
                  log):
    
    # this is a toml file
    config = model.config
    # update the log 
    config["log"] = f"../../0-log/instate_L{0}/instates_level{level}_ST{str(soilthickness).replace('.', '')}.txt"
    
    # casename
    config["casename"] = f"instates_level{level}_ST{str(soilthickness).replace('.', '')}"
    
    # update the reinit
    config["model"]["reinit"] = True
    
    # update the input path_static
    config["input"]["path_static"] = "../staticmaps/staticmaps.nc"
    
    # update the input path_forcing
    config["input"]["path_forcing"] = "../forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc"
    
    # print(config.keys())
    # update the state path_output
    config["state"]["path_output"] = f'instate_level{level}_ST{str(soilthickness).replace(".", "")}.nc'
    config["state"]["path_input"] = None
    
    # print(config.keys())
    
    # Ensure config['input'] is a dictionary
    if not isinstance(config.get("input"), dict):
        config["input"] = {}
        
    # Ensure config['input']['vertical'] is a dictionary
    if not isinstance(config["input"].get("vertical"), dict):
        config["input"]["vertical"] = {}
        
    # Ensure config['input']['vertical']['soilthickness'] is a dictionary
    if not isinstance(config["input"]["vertical"].get("soilthickness"), dict):
        config["input"]["vertical"]["soilthickness"] = {}
        
    # Ensure config['input']['vertical']['soilthickness']['netcdf'] is a dictionary
    if not isinstance(config["input"]["vertical"]["soilthickness"].get("netcdf"), dict):
        config["input"]["vertical"]["soilthickness"]["netcdf"] = {}
        
    # Ensure config['input']['vertical']['soilthickness']['netcdf']['variable'] is a dictionary
    if not isinstance(config["input"]["vertical"]["soilthickness"]["netcdf"].get("variable"), dict):
        config["input"]["vertical"]["soilthickness"]["netcdf"]["variable"] = {}
        
    # Update the input vertical soilthickness
    config["input"]["vertical"]["soilthickness"]["netcdf"]["variable"]["name"] = ST_key
    config["input"]["vertical"]["soilthickness"]["scale"] = soilthickness
    
    # Ensure config['input']['vertical']['soilminthickness'] is a dictionary
    if not isinstance(config["input"]["vertical"].get("soilminthickness"), dict):
        config["input"]["vertical"]["soilminthickness"] = {}
        
    # Ensure config['input']['vertical']['soilminthickness']['netcdf'] is a dictionary
    if not isinstance(config["input"]["vertical"]["soilminthickness"].get("netcdf"), dict):
        config["input"]["vertical"]["soilminthickness"]["netcdf"] = {}
        
    # Ensure config['input']['vertical']['soilminthickness']['netcdf']['variable'] is a dictionary
    if not isinstance(config["input"]["vertical"]["soilminthickness"]["netcdf"].get("variable"), dict):
        config["input"]["vertical"]["soilminthickness"]["netcdf"]["variable"] = {}
    
    # Update the input vertical soilminthickness
    config["input"]["vertical"]["soilminthickness"]["netcdf"]["variable"]["name"] = ST_key
    config["input"]["vertical"]["soilminthickness"]["scale"] = soilthickness
    
    # Update the starttime
    config["starttime"] = start
    
    # Update the endtime
    config["endtime"] = end
    
    #we dont need the output from this period, csv outputs of generic name compete 
    #stopping them from running in parallel
    config.pop("csv", None)
    # print(config.keys())
    
    # Write the config to a TOML file
    output_path = Path(root) / "instates" / f"wflow_sbm_getinstate_level{level}_ST{str(soilthickness).replace('.', '')}.toml"
    with open(output_path, 'w', encoding='utf-8') as toml_file:
        toml.dump(config, toml_file)

    l.info(f"Writing instate file for soil thickness: {soilthickness} to {output_path}")
    return model



if __name__ == "__main__":
    
    l = setup_logging('data/0-log', '01-create_instate_tomls.log')
    try:
        config_fn = snakemake.input.config_fn
        root= snakemake.params.root
        ST_values = list(snakemake.params.ST_values)
        ST_key = snakemake.params.ST_key
        level = snakemake.params.level
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
        for st in ST_values:
            mod = read_model(config_fn, l)
            l.info(f'original mod.config.input: {mod.config["input"]}')
            change_config(model=mod,
                        root=root, 
                        level=level, 
                        ST_key=ST_key,
                        soilthickness=st, 
                        start=start, 
                        end=end,
                        log=l)
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e
    
        
    
        
