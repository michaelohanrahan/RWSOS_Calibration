#%% TODO: 1) how to access best_10params: list[Path | str]. 2) how to save the random select paramset results

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

def select_random_paramset(df, gauge):
    """
    df: csv that contains the best_10params
    gauge: the gauge id
    """
    row = df[df.index == gauge]
    if row.empty:
        raise ValueError("Gauge not found in DataFrame")
    random_column = random.choice(row.columns)
    random_paramset = row[random_column].apply(eval).values[0]

    return random_column, random_paramset


def main(
    l,
    p: Path | str,
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    best_10params: list[Path | str], #TODO: a list of best_10params.csv from previous levels
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
    if not os.path.exists(best_10params) or level == 'level0':
        l.info(f"Best params file not found or level is level0, skipping upper level random params")
    else:
        upgauge_ids = graph[level]["deps"]
        upgauge_int = [int(item) for item in upgauge_ids]
        l.info(f"Updating the following upstream gauges: {upgauge_int}")
        random_sel_results = []
        # for each upstream gauge, set random param from best_10params.csv
        for upgauge in upgauge_int:
            vds_upgauge = vds[vds.value == upgauge]
            mask_up = ds.raster.geometry_mask(vds_upgauge)
            
            #TODO: find the csv that contain this upgauge from best_10params => lets call it csv for now
            # select one from top1 to top10 for this upgauge
            random_column, random_paramset = select_random_paramset(csv, upgauge)
            random_sel_results.append({'upgauges': upgauge, 'selected_column': random_column, 'selected_paramset': random_paramset})
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
                
        df_random_sel_results = pd.DataFrame(random_sel_results)  #TODO: how to save this?
    
    ds.to_netcdf(out)
    l.info(f"Writing dataset to {out}")

    #make sure the lakes h-q rating relation is also copied
    for lin, lout in zip(lakes_in, lakes_out):
        if not os.path.exists(lout):
            shutil.copy(lin, lout)
            l.info(f"Copying {lin} to {lout}")


if __name__ == "__main__":
    
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA')
    csv = pd.read_csv(work_dir / 'best_10params.csv', index_col='gauges')
    csv.index=csv.index.astype(int)
    
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

