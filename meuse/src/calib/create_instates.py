import snakemake
from pathlib import Path
from hydromt_wflow import WflowModel
from hydromt.log import setuplog


def read_model(config_fn):
    # read the model configuration
    logger = setuplog("build", log_level=20)
    model = WflowModel(root=root, mode="r", config_fn=config_fn,  logger=logger)
    model.read_config()
    model.read_grid()
    return model

def change_config(model, 
                  root,
                  level, 
                  ST_key,
                  soilthickness,
                  start,
                  end):
    
    # change the soil thickness in the model configuration
    model.set_config("reinit", False)
    
    #declare it as a nc variable in the 
    model.set_config('input', 'vertical', 'soilthickness', 'netcdf', 'variable', {'name': ST_key})
    model.set_config('input', 'vertical', 'soilminthickness', 'netcdf', 'variable', {'name': ST_key}) 
       
    model.set_config('input','vertical','soilthickness',{"scale":soilthickness})
    model.set_config('input','vertical','soilminthickness',{"scale":soilthickness})
    
    model.set_config("starttime", start)
    model.set_config("endtime", end)     
    
    model.set_config("state", {"path_output":f'instate_level{level}_ST{str(soilthickness).replace(".", "")}.nc'})
    model.set_root(Path(root, "instates"), mode="w+")
    
    model.write_config(f"wflow_sbm_getinstate_level{level}_ST{str(soilthickness).replace('.', '')}.toml")
    return model

if __name__ == "__main__":
    try:
        config_fn = snakemake.input.config_fn
        root= snakemake.params.root
        ST_values = snakemake.params.ST_values
        ST_key = snakemake.params.ST_key
        level = snakemake.params.level
        start = snakemake.params.starttime
        end = snakemake.params.endtime
        os.makedirs("instates", exist_ok=True)
    except:
        os.chdir(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse')
        print("Running script in test mode")
        config_fn= Path("wflow_sbm.toml")   
        root= Path("data/3-input")
        ST_values =[999, 9999]
        ST_key = "wrongsoilthicknesskey"
        level = 0
        start = "2005-01-01T00:00:00"
        end = "2005-03-01T00:00:00"
        os.makedirs("instates", exist_ok=True)
    
    
    for st in ST_values:
        mod = read_model(config_fn)
        change_config(model=mod,
                      root=root, 
                      level=level, 
                      ST_key=ST_key,
                      soilthickness=st, 
                      start=start, 
                      end=end)
        
    
        
    
        
