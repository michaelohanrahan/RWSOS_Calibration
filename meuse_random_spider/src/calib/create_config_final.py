# import snakemake
from pathlib import Path
from hydromt_wflow import WflowModel
from hydromt.log import setuplog
import traceback
from setuplog import setup_logging
import os
import toml

def add_nc_variable(wflow_string, var_name, out_cfg):
    """
    Add a nested structure for netCDF variables in the TOML config.
    """
    keys = wflow_string.split('.')
    current = out_cfg
    for key in keys[:-1]:
        if key not in current:
            current[key] = {}
        current = current[key]
    
    if keys[-1] not in current:
        current[keys[-1]] = {"name": var_name}
    return out_cfg


def add_scale_offset(nc_wflow_string, out_cfg, scale=1.0, offset=0.0):
    """
    Add scale and offset to a nested structure in the TOML config.
    """
    keys = nc_wflow_string.split('.')
    current = out_cfg
    for key in keys:
        if key not in current:
            current[key] = {}
        current = current[key]
    
    current["scale"] = scale
    current["offset"] = offset
    return out_cfg


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
                  gaugemap,
                  l,
                  outcfg):
    # this is a toml file
    config = model.config
    l.info(F"setting up final model run config for {topx}")
    # update the log 
    config["log"] = f"../../0-log/final_model_run_{topx}/final_model_run_{topx}.txt"
    
    # casename
    config["casename"] = f"final_model_run_{topx}"
    
    # update the reinit
    config["model"]["reinit"] = False
    
    # update the input path_static
    config["input"]["path_static"] = f"staticmaps.nc"
    
    # update the input path_forcing
    config["input"]["path_forcing"] = forcing_path.as_posix()
    
    # update the state path_output
    config["state"]["path_input"] = f"instates/instate_final.nc"
    config["state"]["path_output"] = f'outstate_{topx}.nc'

    # Update the starttime
    config["starttime"] = start
    
    # Update the endtime
    config["endtime"] = end
    
    #we dont need the output from this period, csv outputs of generic name compete 
    #stopping them from running in parallel
    if "csv" in config:
        config.pop("csv", None)
        
    if "netcdf" not in config:
        config["netcdf"] = {}
    
    config["netcdf"]["path"] = f"../../4-output/output_{topx}/output_scalar.nc"
    config["netcdf"]["variable"] = [
        {
            "name":"Q",
            "map":gaugemap,
            "parameter":"lateral.river.q_av"
        },
        {
            "name":"Q_hbv",
            "map":"gauges_hbv",
            "parameter":"lateral.river.q_av"
        }
    ]
    
    # look for the input.lateral.river.floodplain.n.netcdf.variable
    # ... if not add it "N_Floodplain"
    config = add_nc_variable("input.lateral.river.flooplain.n.netcdf.variable", 
                              "N_Floodplain",
                              out_cfg=config)
    
    config = add_scale_offset("input.lateral.river.flooplain.n", 
                               scale=1.0,
                               offset=0.0,
                               out_cfg=config)
    
    # Write the config to a TOML file
    output_path = outcfg
    with open(output_path, 'w', encoding='utf-8') as toml_file:
        toml.dump(config, toml_file)
    l.info(f"Updated config saved to {output_path}")
    return model



if __name__ == "__main__":
    
    l = setup_logging('data/0-log', f'02-final_tomls_{snakemake.params.topx}.log')
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        wflow_mod = read_model(mod.params.root, mod.input.config_fn)
        change_config(model=wflow_mod,
                    start=mod.params.starttime, 
                    end=mod.params.endtime,
                    forcing_path=mod.params.forcing_path,
                    topx=mod.params.topx,
                    gaugemap=mod.params.gaugemap,
                    l=l,
                    outcfg=mod.output.cfg)
    
    
    
    # try:
    #     config_fn = snakemake.input.config_fn
    #     root= snakemake.params.root
    #     level = snakemake.params.level
    #     start = snakemake.params.starttime
    #     end = snakemake.params.endtime
    #     os.makedirs("instates", exist_ok=True)
        
    
    # except Exception as e:
    #     l.error(f"An error occurred importing params: {e}")
    #     l.error(traceback.format_exc())
    #     raise e

    # try:
    #     mod = read_model(config_fn, l)
    #     change_config(model=mod,
    #                 root=root, 
    #                 level=level, 
    #                 start=start, 
    #                 end=end,
    #                 log=l,
    #                 outcfg=snakemake.output.cfg)
    # except Exception as e:
    #     l.error(f"An error occurred: {e}")
    #     l.error(traceback.format_exc())
    #     raise e
    
        
    
        
