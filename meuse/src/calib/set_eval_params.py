from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from hydromt.raster import RasterDataset


def main(
    params: Path | str,
    staticmaps: Path| str,
    sub_catch: Path | str,
    params_lname: tuple | list,
    params_method: tuple | list,
    out: Path | str,
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
    
    par_da = [
        ds[var] for var in params_lname
    ]

    for gauge in params_ds.index:
        select_vds = vds[vds.value.isin([gauge])]
        mask = ds.raster.geometry_mask(select_vds)
        for idx, value in enumerate(params_ds.loc[gauge,:].values):
            if params_method[idx] == "mult":
                par_da[idx].values[mask] *= value
            elif params_method[idx] == "set":
                par_da[idx].values[mask] = value

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