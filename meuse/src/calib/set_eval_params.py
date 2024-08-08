from pathlib import Path
from setuplog import setup_logging
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydromt.raster import RasterDataset
import random
import ast



def main(
    params: Path | str,
    staticmaps: Path| str,
    sub_catch: Path | str,
    params_lname: tuple | list,
    params_method: tuple | list,
    level: int,
    out: Path | str, #nc
    txt: Path | str, #txt
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
    prev_level = level-1
    prev_level=f"L{prev_level}"
    save_old_dir = out.parent / "intermediate"
    save_old = f"staticmaps_{prev_level}.nc"
    os.makedirs(save_old_dir, exist_ok=True)
    shutil.copy(staticmaps, Path(save_old_dir) / save_old)
    
    """_summary_"""
    # Get the staticmaps
    with xr.open_dataset(staticmaps) as _r:
        ds = _r.load()

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    # Get the relevant sub catchments
    vds = vds.astype({"value": int})

    # Read the parameters
    params_ds = pd.read_csv(params, index_col="gauges")
    params_ds.index = params_ds.index.astype(int)
    
    par_da = [
        ds[var] for var in params_lname
    ]
    
    for gauge, param_set in selected_params.items():
        # Select the sub-catchments corresponding to the current gauge
        select_vds = vds[vds.value == gauge]
        mask = ds.raster.geometry_mask(select_vds)
        
        #param_set
        param_row = params_ds.loc[gauge, "Top_1"]
        param_set = ast.literal_eval(param_row)
        l.info(f"Applying parameters for gauge {gauge}: {param_set}")
        #the

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
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        l = setup_logging('data/0-log', f'08-initial_instate_tomls_L{snakemake.params.level}.log')

        main(
            params = mod.input.best_params,
            staticmaps = mod.params.staticmaps,
            sub_catch = mod.params.sub_catch,
            params_lname = mod.params.params_lname,
            params_method = mod.params.params_method,
            level=mod.params.level,
            out=mod.output.done,
        )

    else:
        main(
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/best_params.csv",
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
            ["KsatHorFrac", "RootingDepth", "SoilThickness"],
            ["set", "mult", "mult"],
            "",
        )
    pass