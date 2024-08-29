from pathlib import Path
from setuplog import setup_logging
from hydromt.raster import RasterDataset
import geopandas as gpd
import pandas as pd
import xarray as xr
import ast
import os
import shutil
import sys

def save_copy(l,level, staticmaps):
    prev_level = level-1
    prev_level=f"L{prev_level}"
    save_old_dir = Path(Path(staticmaps).parent, "intermediate")
    save_old = f"staticmaps_{prev_level}.nc"
    os.makedirs(save_old_dir, exist_ok=True)
    shutil.copy(staticmaps, Path(save_old_dir,save_old))
    assert os.path.exists(Path(save_old_dir,save_old)), f"Failed to copy {staticmaps} to {Path(save_old_dir,save_old)}"
    l.info(f"Copied {staticmaps} to {Path(save_old_dir,save_old)}")
    return Path(save_old_dir,save_old)

def main(
    l,
    params: Path | str,
    staticmaps: Path| str,
    sub_catch: Path | str,
    params_lname: tuple | list,
    params_method: tuple | list,
    level: int,
    out: Path | str, #nc
):
    """
    Apply evaluation parameters to the static maps.

    Args:
        params (Path or str): Path to the best params csv.
        staticmaps (Path or str): Path to the original staticmaps file.
        sub_catch (Path or str): Path to the sub catchments file.
        params_lname (tuple or list): List of parameter names.
        params_method (tuple or list): List of parameter methods.
        out (Path or str): Path to the output file.

    Returns:
        None
    """
    #preseve the precursor gridfiles
    copied = save_copy(l,level, staticmaps)
    assert os.path.exists(copied), f"Failed to copy {staticmaps} to {copied}"
    l.info(f"Copied {staticmaps} to {copied}")
    
    # Get the staticmaps
    with xr.open_dataset(staticmaps) as _r:
        ds = _r.load()

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    
    # Get the relevant sub catchments
    vds = vds.astype({"value": int})

    # Read the parameters
    params_ds = pd.read_csv(params, index_col="gauge")
    params_ds.index = params_ds.index.astype(int)
    
    par_da = [
        ds[var] for var in params_lname
    ]
    
    for gauge in params_ds.index.values:
        # Select the sub-catchments corresponding to the current gauge
        select_vds = vds[vds.value == gauge]
        mask = ds.raster.geometry_mask(select_vds)
        #param_set
        param_row = params_ds.loc[gauge, "Top_1"]
        param_set = ast.literal_eval(param_row)
        l.info(f"Applying the best parameters for gauge {gauge}: {param_set}")
        for idx, param_value in enumerate(param_set.values()):
            if params_method[idx] == "mult":
                par_da[idx].values[mask] *= param_value
            elif params_method[idx] == "set":
                par_da[idx].values[mask] = param_value
            elif params_method[idx] == "add":
                par_da[idx].values[mask] += param_value

    for idx, da in enumerate(par_da):
        ds[params_lname[idx]] = da

    ds.to_netcdf(staticmaps)
    l.info(f"Updated staticmaps saved to {staticmaps}")
    
    with open(out, "w") as _w:
        _w.write("Done!\n")
    pass          


if __name__ == "__main__":
    
    l = setup_logging('data/0-log', f'08-Set_Params.log')

    if "snakemake" in globals():
        mod = globals()["snakemake"]
        main(
            l,
            params = mod.input.best_params,
            staticmaps = mod.params.staticmaps,
            sub_catch = mod.params.sub_catch,
            params_lname = mod.params.params_lname,
            params_method = mod.params.params_method,
            level=mod.params.level,
            out=mod.output.done,
        )

    else:
        if sys.platform == "linux":
            DRIVE = "/p/"
        else:
            DRIVE = "p:/"
        os.chdir(f"{DRIVE}11209265-grade2023/wflow/RWSOS_Calibration/meuse")
        main(
            l,
            Path(f"{DRIVE}11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/best_params.csv"),
            Path(f"{DRIVE}11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/3-input/staticmaps/staticmaps.nc"),
            Path(f"{DRIVE}11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/3-input/staticgeoms/subcatch_Hall.geojson"),
            ["ksathorfrac_BRT_250", "f_", "RootingDepth_obs_15", "SoilThickness_manual_cal", "N_River", "MaxLeakage_manual_cal"],
            ["mult", "mult", "mult", "mult", "mult", "add"],
            0,
            Path(f"{DRIVE}11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/done.txt"),
        )
    pass