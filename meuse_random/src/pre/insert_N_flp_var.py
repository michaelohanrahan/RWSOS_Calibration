import xarray as xr 
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
import argparse
import traceback
import shutil

def main(gridfile_in, config_fn_in, geoms, config_fn_out, gridfile_out, geoms_out):
    # Open the dataset
    ds = xr.open_dataset(gridfile_in)

    # Create a constant array with the same shape and dimensions as 'wflow_dem'
    constant_value = 0.072
    var = np.full_like(ds['wflow_dem'], constant_value)

    # Apply the NaN mask from 'wflow_dem' to the constant array
    var = np.where(np.isnan(ds['wflow_dem']), np.nan, var)

    # Create a DataArray with the same dims, coords, and attrs as 'wflow_dem'
    da = xr.DataArray(var, 
                      dims=ds['wflow_dem'].dims, 
                      coords=ds['wflow_dem'].coords, 
                      attrs=ds['wflow_dem'].attrs)

    # now we insert N_Floodplain into the dataset
    ds = ds.assign(N_Floodplain=da)
    ds['N_Floodplain'].attrs = {'long_name': 'N_Floodplain', 
                                'units': '-'}
    ds.to_netcdf(gridfile_out)

    # copy the geoms
    shutil.copytree(geoms, geoms_out)
    shutil.copy(config_fn_in, config_fn_out)
    

if __name__ == "__main__":
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            gridfile_in = mod.input.gridfile_in
            config_fn_in = mod.input.config_fn_in
            geoms = mod.input.geoms
            config_fn_out = mod.output.config_fn_out
            gridfile_out = mod.output.gridfile_out
            geoms_out = mod.output.geoms_out
        else:
            raise ValueError("Snakemake not in global scope")
            # gridfile_in = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps_old.nc"
            # config_fn_in = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\addksathorfrac\wflow_sbm_addksathorfrac.toml"
            # geoms = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\addksathorfrac\staticgeoms"
            # config_fn_out = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\addflp\wflow_sbm_addflp.toml"
            # gridfile_out = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\addflp\staticmaps\staticmaps.nc"
            # geoms_out = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\addflp\staticgeoms"

        main(gridfile_in, config_fn_in, geoms, config_fn_out, gridfile_out, geoms_out)
    except Exception as e:
        traceback.print_exc()
        print(e)
        raise e
        
        
=======

# Open the dataset
ds = xr.open_dataset(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps_old.nc")

# Create a constant array with the same shape and dimensions as 'wflow_dem'
constant_value = 0.072
var = np.full_like(ds['wflow_dem'], constant_value)

# Apply the NaN mask from 'wflow_dem' to the constant array
var = np.where(np.isnan(ds['wflow_dem']), np.nan, var)

# Create a DataArray with the same dims, coords, and attrs as 'wflow_dem'
da = xr.DataArray(var, dims=ds['wflow_dem'].dims, coords=ds['wflow_dem'].coords, attrs=ds['wflow_dem'].attrs)

# now we insert N_Floodplain into the dataset
ds = ds.assign(N_Floodplain=da)
ds['N_Floodplain'].attrs = {'long_name': 'N_Floodplain', 'units': '-'}
ds.to_netcdf(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps.nc")
>>>>>>> 50f0b2a433687946660e04c31a5fc4d5b3a8b94a
