from pathlib import Path
from setuplog import setup_logging
import shutil
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydromt.raster import RasterDataset
import random


def main(
    params: Path | str,
    staticmaps: Path| str,
    sub_catch: Path | str,
    params_lname: tuple | list,
    params_method: tuple | list,
    level: int,
    out: Path | str,
):
    """
    Apply evaluation parameters to the static maps.

    Args:
        params (Path or str): Path to the best 10 params csv. (e.g., "best_10params.csv")
        staticmaps (Path or str): Path to the original staticmaps file. #TODO: there will be 10 original staticmaps
        sub_catch (Path or str): Path to the sub catchments file.
        params_lname (tuple or list): List of parameter names.
        params_method (tuple or list): List of parameter methods.
        out (Path or str): Path to the output file.

    Returns:
        None
    """
    #preseve the precursor gridfiles
    #TODO: there will be 10 original staticmaps
    prev_level = level-1
    if prev_level < 0:
        prev_level = 'original'
    else:
        f"L{prev_level}"
    save_old_dir = out.parent / "intermediate"
    save_old = f"staticmaps_{prev_level}.nc"
    
    shutil.copy(staticmaps, Path(save_old_dir) / save_old)  # save an old staticmaps
    
    # """_summary_"""
    # # Get the staticmaps
    # with xr.open_dataset(staticmaps) as _r:
    #     ds = _r.load()

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    # Get the relevant sub catchments
    vds = vds.astype({"value": int})

    # Read the parameters
    # Index by gauges
    params_ds = pd.read_csv(params, index_col="gauges")
    params_ds.index = params_ds.index.astype(int)
    
    #TODO: modify for later levels where there already multiple original staticmaps
    for column in params_ds.columns:
        # Get the staticmaps for each parameter set, TODO: select the right staticmaps (the right staticmaps needs to be registered in set_calib_params.py)
        with xr.open_dataset(staticmaps) as _r:
            ds = _r.load()

        selected_params = params_ds[column].apply(eval)  # Convert string representations of dicts to dicts

        par_da = [ds[var] for var in params_lname]

        for gauge, param_set in selected_params.items():
            # Select the sub-catchments corresponding to the current gauge
            select_vds = vds[vds.value == gauge]
            mask = ds.raster.geometry_mask(select_vds)

            for idx, param_value in enumerate(param_set.values()):
                if params_method[idx] == "mult":
                    par_da[idx].values[mask] *= param_value
                elif params_method[idx] == "set":
                    par_da[idx].values[mask] = param_value
                elif params_method[idx] == "add":
                    par_da[idx].values[mask] += param_value

        for idx, da in enumerate(par_da):
            ds[params_lname[idx]] = da

        # Save each modified staticmap with a unique name
        modified_staticmap_path = out / f"staticmaps_{column}.nc"
        ds.to_netcdf(modified_staticmap_path)
    
    with open(out, "w") as _w:
        _w.write("Done!\n")
    pass          


if __name__ == "__main__":
    
    work_dir = Path(r"c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA")
    
    import pickle as pk
    with open(work_dir/'create_set_params.pkl', 'rb') as f:
        dict = pk.load(f)
    lnames = dict['lnames']
    methods = dict['methods']
    df = dict['ds']
    
    # inputs
    params = work_dir / "best_10params.csv"
    staticmaps = work_dir / "staticmaps/staticmaps.nc"
    sub_catch = work_dir / 'subcatch_Hall.geojson'
    params_lname = lnames
    params_method = methods
    level = 10
    out = work_dir / "staticmaps"
    
    main(params, staticmaps, sub_catch, params_lname, params_method, level, out)
    
    
    # if "snakemake" in globals():
    #     mod = globals()["snakemake"]
    #     l = setup_logging('data/0-log', f'08-initial_instate_tomls_L{snakemake.params.level}.log')
    #     main(
    #         mod.input.best_params,
    #         mod.params.staticmaps,
    #         mod.params.sub_catch,
    #         mod.params.params_lname,
    #         mod.params.params_method,
    #         mod.params.level,
    #         mod.output.done,
    #     )

    # else:
    #     main(
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/best_params.csv",
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
    #         ["KsatHorFrac", "RootingDepth", "SoilThickness"],
    #         ["set", "mult", "mult"],
    #         "",
    #     )
    # pass