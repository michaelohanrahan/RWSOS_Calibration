from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydromt.raster import RasterDataset
import random

#TODO: to be tested

def main(
    params: Path | str,
    staticmaps: Path| str,
    sub_catch: Path | str,
    params_lname: tuple | list,
    params_method: tuple | list,
    out: Path | str,
):
    """_summary_"""
    # Get the staticmaps
    with xr.open_dataset(staticmaps) as _r:
        ds = _r.load()

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    # Get the relevant sub catchments
    vds = vds.astype({"value": int})

    # Read the parameters
    # Index by gauges
    params_ds = pd.read_csv(params, index_col="gauges")
    params_ds.index = params_ds.index.astype(int)
    
    # Randomly select one of the Top_x columns
    selected_column = random.choice(params_ds.columns)
    selected_params = params_ds[selected_column].apply(eval)  # Convert string representations of dicts to dicts
    
    par_da = [
        ds[var] for var in params_lname
    ]
    
    
    #TODO: params_method: "add"; co-scaling n_land, n_
    for gauge, param_set in selected_params.items():
        # Select the sub-catchments corresponding to the current gauge
        select_vds = vds[vds.value == gauge]
        mask = ds.raster.geometry_mask(select_vds)

        for idx, param_value in enumerate(param_set.values()):
            if params_method[idx] == "mult":
                par_da[idx].values[mask] *= param_value
            elif params_method[idx] == "set":
                par_da[idx].values[mask] = param_value
    

    for idx, da in enumerate(par_da):
        ds[params_lname[idx]] = da

    ds.to_netcdf(staticmaps)

    with open(out, "w") as _w:
        _w.write("Done!\n")
    pass          


if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        main(
            mod.input.best_params,
            mod.params.staticmaps,
            mod.params.sub_catch,
            mod.params.params_lname,
            mod.params.params_method,
            mod.output.done,
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