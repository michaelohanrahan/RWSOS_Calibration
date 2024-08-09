#%%TODO: a path to random_df created in previous level

#%%
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
    p: Path | str,
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    random_df: Path | str, #TODO: a path to random_df created in previous rule
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
    # Load original staticmaps
    ds = xr.open_dataset(p)
    
    # Load the random_df
    random_df = pd.read_csv(random_df, index_col=0)

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    vds = vds.astype({"value": int})
    
    # Set param for current level gauges
    gauge_ids = graph[level]["elements"]
    gauge_int = [int(item) for item in gauge_ids]
    l.info(f"Updating the following current level gauges: {gauge_int}")
    vds_current = vds[vds.value.isin(gauge_int)]
    mask = ds.raster.geometry_mask(vds_current)
    for idx, (key, value) in enumerate(params.items()):
        vds[key] = value   # what is this line used for?
        da = ds[params_lname[idx]]
        if params_method[idx] == "mult":
            da.values[mask] *= value
        elif params_method[idx] == "set":
            da.values[mask] = value
        elif params_method[idx] == "add":
            da.values[mask] += value
        ds[params_lname[idx]] = da
    
    # Set param for upstream gauges
    if not os.path.exists(random_df) or level == 'level0':
        l.info(f"random_df file not found or level is level0, skipping upper level random params")
    else:
        # select all the relevant upstream gauges (including further upstream)
        upgauge_ids = graph[level]['deps']
        upgauge_int = [int(item) for item in upgauge_ids]
        l.info(f"Updating the following upstream gauges: {upgauge_int}")
        
        # select matching row from random_df
        _mask = pd.Series([True] * len(random_df))
        for key, value in params.items():
            _mask = _mask & (random_df[key] == value)
        random_df_sel = random_df[_mask]
        # check if the matching is unique
        if len(random_df_sel) != 1:
            raise ValueError(f"Matching row in random_df is not unique! {len(random_df_sel)} rows found.")
        
        # for each upstream gauge, set random param from random_df_sel
        for upgauge in upgauge_int:
            vds_upgauge = vds[vds.value == upgauge]
            mask_up = ds.raster.geometry_mask(vds_upgauge)
            # select random paramset
            random_paramset = random_df_sel[str(float(upgauge))].apply(eval)[0]
            # modify staticmap based on the random paramset
            for idx, value in enumerate(random_paramset.values()):
                da = ds[params_lname[idx]]
                if params_method[idx] == "mult":
                    da.values[mask_up] *= value
                elif params_method[idx] == "set":
                    da.values[mask_up] = value
                elif params_method[idx] == "add":
                    da.values[mask_up] += value          
                ds[params_lname[idx]] = da
    
    ds.to_netcdf(out)
    l.info(f"Writing dataset to {out}")

    #make sure the lakes h-q rating relation is also copied
    for lin, lout in zip(lakes_in, lakes_out):
        if not os.path.exists(lout):
            shutil.copy(lin, lout)
            l.info(f"Copying {lin} to {lout}")


if __name__ == "__main__":
    
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA')
    random_df = pd.read_csv(work_dir / 'random_df.csv', index_col=0)
    
    p = work_dir / 'staticmaps/staticmaps.nc'
    
    import pickle as pk
    with open(work_dir/'create_set_params.pkl', 'rb') as f:
        dict = pk.load(f)
    params_lname = dict['lnames']
    params_method = dict['methods']
    df = dict['ds']
    
    from snakemake.utils import Paramspace
    paramspace = Paramspace(df)
    params = paramspace.instance
    
    level = 'level5'
    import json
    graph = json.load(open(Path(work_dir / 'Hall_levels_graph.json')))
    
    sub_catch = work_dir / 'subcatch_Hall.geojson'
    vds = gpd.read_file(sub_catch)
    vds = vds.astype({"value": int})
    
    gauge_ids = graph['level4']["elements"]
    upgauge_new = graph['level4']['deps']
    


    
    
    # l=setup_logging('data/0-log', '03-set_calib_params.log')
    # try:
    #     if "snakemake" in globals():
    #         mod = globals()["snakemake"]

    #         main(
    #             l,
    #             mod.params.best_params,
    #             mod.params.dataset,
    #             mod.params.params,
    #             mod.params.params_lname,
    #             mod.params.params_method,
    #             mod.params.level,
    #             mod.params.graph[mod.params.level]["elements"],
    #             mod.params.sub_catch,
    #             mod.params.lake_in,
    #             mod.output.lake_hqs, 
    #             mod.output.staticmaps,
    #         )

    #     else:
    #         main(
    #             "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
    #             {"ksat": 5, "rd": 0.5, "st": 0.5},
    #             ["KsatHorFrac", "RootingDepth", "SoilThickness"],
    #             ["set", "mult", "mult"],
    #             ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000'],
    #             "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
    #             "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/out.nc"
    #         )
    #     pass

    # except Exception as e:
    #     l.error(f"An error occurred: {e}")
    #     l.error(traceback.format_exc())
    #     raise e

