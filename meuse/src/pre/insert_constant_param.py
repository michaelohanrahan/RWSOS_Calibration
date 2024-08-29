
import xarray as xr 
import numpy as np
from pathlib import Path
import os
import sys
sys.path.append(os.getcwd())
from src.calib.setuplog import setup_logging
import traceback 

def main(logger, 
         staticmaps, 
         constant_value:float=None,
         units:str=None,
         long_name:str='N_Floodplain',  
         out:str|Path='data/3-input/staticmaps.nc',
         ):
    logger.info(f"cwd: {os.getcwd()}")
    cwd = Path(os.getcwd()).as_posix()
    logger.info(f"Adding {long_name} to staticmaps.nc")
    logger.info(f"Copying from {staticmaps} to {out}")
    out = Path(out)
    c_dict = {
            'N_Floodplain':
                {
                'value':0.072,
                'units':'-'
                },
        }
    if constant_value is None:
        if long_name in c_dict.keys():
            constant_value = c_dict[long_name]['value']
            units = c_dict[long_name]['units']
        else:
            logger.error(f"Constant value for {long_name} not found in the dictionary. Please provide a value.")
            raise ValueError(f"Constant value for {long_name} not found in the dictionary. Please provide a value.")
    
    if units is None:
        logger.error(f"Units for {long_name} not provided. Please provide units.")
        raise ValueError(f"Units for {long_name} not provided. Please provide units.")
    
    # Open the dataset
    ds = xr.open_dataset(Path(cwd,staticmaps))

    # Create a constant array with the same shape as 'wflow_dem'
    var = np.full_like(ds['wflow_dem'], constant_value)

    # Apply the NaN mask from 'wflow_dem' to the constant array
    var = np.where(np.isnan(ds['wflow_dem']), np.nan, var)

    # Create a DataArray with the same dims, coords, and attrs as 'wflow_dem'
    da = xr.DataArray(var, dims=ds['wflow_dem'].dims, coords=ds['wflow_dem'].coords, attrs=ds['wflow_dem'].attrs)

    # Insert N_Floodplain into the dataset
    ds = ds.assign(N_Floodplain=da)
    ds[long_name].attrs = {'long_name': long_name, 
                                'units': units}

    # Save the updated dataset
    ds.to_netcdf(Path(cwd, out))

    # Log the completion
    logger.info(f"Updated staticmaps.nc with N_Floodplain. Output saved to {out}")
    # Indicate that processing is done by writing to the output file
    Path(Path(cwd, out).parent, "added_vars.txt").write_text(f"added {long_name} to staticmaps.nc")

if __name__ == "__main__":
    
    l = setup_logging('data/0-log', '08-Set_Params.log')
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            main(
                l,
                staticmaps = mod.params.staticmaps,
                params_lname = mod.params.params_lname,
                params_method = mod.params.params_method,
                out=mod.output.done,
            )
        else:
            os.chdir('/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse')
            current_dir = Path(os.getcwd())
            l.info(f"{'*'*10} Running script from {current_dir} {'*'*10}")
            if str(current_dir.name) not in ['meuse', 'rhine']:
                l.error(f"Please run this script from the 'meuse' or 'rhine' directory, not {str(current_dir.name)}")
                raise ValueError(f"Please run this script from the 'meuse' or 'rhine' directory, not {os.getcwd()}")
            main(
                l,
                staticmaps = Path("data/2-interim/addksathorfrac/staticmaps/staticmaps.nc"),
                long_name = 'N_Floodplain',
            )
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e