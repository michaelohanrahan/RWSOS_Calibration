from pathlib import Path
import pandas as pd
import geopandas as gpd
import xarray as xr
from hydromt.raster import RasterDataset
from setuplog import setup_logging
import shutil
import traceback
import os
import random
import numpy as np


def main(
    l,
    p: list[Path | str],  # a list of 10 staticmaps paths from the previous level
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    level: str,
    graph: dict,
    sub_catch: Path | str,
    lakes_in: Path | str,
    lakes_out: Path | str,
    out: Path | str,
):
    """
    Load the dataset (staticmaps) then update params for the desired gauge_ids
    the paramslname should be the same as the dataset key
    """
    gauge_ids = graph[level]["elements"]
    
    # # Load the current main dataset
    # ds = xr.open_dataset(p)
    # l.info(f"Loading dataset from {p}")
    
    # Make integers of the gauge_ids
    gauge_int = [
        int(item) for item in gauge_ids
    ]
    l.info(f"Updating the following gauge_ids: {gauge_int}")

    # Load the geometries
    # vds
    vds = gpd.read_file(sub_catch)
    
    # Get the relevant sub catchments
    #vds_d is now a geodataframe
    vds_d = vds_d.astype({"value": int})
    #vds_d is now a geodataframe with only the gauges of interest
    vds_d = vds_d[vds_d.value.isin(gauge_int)]
    
    # # Create the data mask
    # mask = ds.raster.geometry_mask(vds_d)

    # Loop through the parameters
    for idx, (key, value) in enumerate(params.items()):
        # randomly select a staticmap
        selected_staticmap_path = random.choice(p)
        ds = xr.open_dataset(selected_staticmap_path)
        l.info(f"Loading dataset from {selected_staticmap_path}")
        
        # Create the mask based on the geometries
        mask = ds.raster.geometry_mask(vds_d)
        
        #The key is the var id
        vds_d[key] = value
        da = ds[params_lname[idx]]
        if params_method[idx] == "mult":
            da.values[mask] *= value
        elif params_method[idx] == "set":
            da.values[mask] = value
        elif params_method[idx] == "add":
            da.values[mask] += value          
        ds[params_lname[idx]] = da
        
        # Save the modified staticmap, TODO: save one or list of all modified staticmaps?
        ds.to_netcdf(out)
        l.info(f"Writing dataset to {out}")

    #make sure the lakes h-q rating relation is also copied
    for lin, lout in zip(lakes_in, lakes_out):
        if not os.path.exists(lout):
            shutil.copy(lin, lout)
            l.info(f"Copying {lin} to {lout}")


if __name__ == "__main__":
    
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA')
    
    import pickle as pk
    with open(work_dir/'create_set_params.pkl', 'rb') as f:
        dict = pk.load(f)
    lnames = dict['lnames']
    methods = dict['methods']
    df = dict['ds']
    
    from snakemake.utils import Paramspace
    paramspace = Paramspace(df)
    params = paramspace.instance
    
    
    
    
    
    
    l=setup_logging('data/0-log', '03-set_calib_params.log')
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]

            main(
                l,
                mod.params.best_params,
                mod.params.dataset,
                mod.params.params,
                mod.params.params_lname,
                mod.params.params_method,
                mod.params.level,
                mod.params.graph[mod.params.level]["elements"],
                mod.params.sub_catch,
                mod.params.lake_in,
                mod.output.lake_hqs, 
                mod.output.staticmaps,
            )

        else:
            main(
                "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
                {"ksat": 5, "rd": 0.5, "st": 0.5},
                ["KsatHorFrac", "RootingDepth", "SoilThickness"],
                ["set", "mult", "mult"],
                ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000'],
                "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
                "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/out.nc"
            )
        pass

    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e

