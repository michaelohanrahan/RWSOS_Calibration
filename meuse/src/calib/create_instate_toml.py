import snakemake
from pathlib import Path
from hydromt_wflow import WflowModel
from hydromt.log import setuplog
import traceback
from setuplog import setup_logging

def read_model(config_fn, log):
    # read the model configuration
    model = WflowModel(root=root, mode="r", config_fn=config_fn,  logger=log)
    model.read_config()
    model.read_grid()
    return model

def change_config(model, 
                  root,
                  level, 
                  ST_key,
                  soilthickness,
                  start,
                  end,
                  log):
    
    # change the soil thickness in the model configuration
    model.set_config("reinit", False)
    
    #declare it as a nc variable in the 
    model.set_config('input', 'vertical', 'soilthickness', 'netcdf', 'variable', {'name': ST_key})
    model.set_config('input', 'vertical', 'soilminthickness', 'netcdf', 'variable', {'name': ST_key}) 
       
    model.set_config('input','vertical','soilthickness',{"scale":soilthickness})
    model.set_config('input','vertical','soilminthickness',{"scale":soilthickness})
    
    model.set_config('input', {'path_forcing':"../forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc"})
    model.set_config('input', {'path_static':"../staticmaps/staticmaps.nc"})
    model.set_config("starttime", start)
    model.set_config("endtime", end)     
    
    model.set_config("state", {"path_output":f'instate_level{level}_ST{str(soilthickness).replace(".", "")}.nc'})
    model.set_root(Path(root, "instates"), mode="w+")
    
    l.info(f"Writing instate file for soil thickness: {soilthickness}")
    model.write_config(f"wflow_sbm_getinstate_level{level}_ST{str(soilthickness).replace('.', '')}.toml")
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
        l.info(f"Creating instate files for soil thickness values: {ST_values}")
    except:
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
        
    
    
    try:
        for st in ST_values:
            mod = read_model(config_fn, l)
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
    
        
    
        
