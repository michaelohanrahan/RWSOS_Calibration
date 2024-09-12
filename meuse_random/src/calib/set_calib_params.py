
from pathlib import Path
import pandas as pd
import geopandas as gpd
import xarray as xr
from hydromt.raster import RasterDataset
from setuplog import setup_logging
import shutil
import traceback
import os
from latin_hyper_paramspace import create_set_all_levels
import random
import numpy as np
import json
from icecream import ic 
import ast

def handle_cs(params_sname, params_lname, params_method):
    '''
    Handle comma-separated items in params_sname and params_lname,
    adjusting params_method accordingly.
    '''
    new_sname = []
    new_lname = []
    new_method = []

    for sname, lname, method in zip(params_sname, params_lname, params_method):
        # #ic(sname, lname, method)
        if ',' in lname:
            lnames = lname.split(',')
            snames = sname.split(',')
            new_lname.extend(lnames)
            new_sname.extend(snames)
            new_method.extend([method] * len(lnames))
        else:
            new_lname.append(lname)
            new_sname.append(sname)
            new_method.append(method)

    return new_sname, new_lname, new_method
    
def main(
    l,
    p: Path | str,
    params: dict,
    params_sname: tuple | list,
    params_lname: tuple | list,
    params_method: tuple | list,
    random_params: Path | str, 
    level: str,
    graph: dict,
    sub_catch: Path | str,
    lakes_in: Path | str,
    lakes_out: Path | str,
    out: Path | str,
):
    """
    Load the dataset (staticmaps) then update params for the desired gauge_ids
    the params lname should be the same as the dataset key
    """
    #handle comma separated params
    params_sname, params_lname, params_method = handle_cs(list(params_sname), list(params_lname), list(params_method))
    #dict of sname to method
    params_sname_to_method = {sname: method for sname, method in zip(params_sname, params_method)}
    #dict of sname to lname
    params_sname_to_lname = {sname: lname for sname, lname in zip(params_sname, params_lname)}
    #ic(params_sname_to_method)
    #ic(params_sname_to_lname)
    # Load original staticmaps
    ds = xr.open_dataset(p)
    
    # Load the geometries
    vds = gpd.read_file(sub_catch)
    vds = vds.astype({"value": int})
    
    # Set param for current level gauges
    gauge_ids = graph[level]["elements"]
    gauge_int = [int(item) for item in gauge_ids]
    # l.info(f"Updating the following current level gauges: {gauge_int}")
    
    vds_current = vds[vds.value.isin(gauge_int)]
    mask = ds.raster.geometry_mask(vds_current)
    
    # Go through all params and set the gauges
    for key, value in params.items():
        da = ds[params_sname_to_lname[key]]
        
        if params_sname_to_method[key] == "mult":
            da.values[mask] *= value
            
        elif params_sname_to_method[key] == "set":
            da.values[mask] = value
            
        elif params_sname_to_method[key] == "add":
            da.values[mask] += value
            
        ds[params_sname_to_lname[key]] = da
    
    # Set param for upstream gauges
    if level == 'level0':
        l.info(f"random_params file not found or level is level0, skipping upper level random params")

    else:
         # Load the random_params
        random_params = pd.read_csv(random_params, index_col=0) 
        # select all the relevant upstream gauges (including further upstream)
        upgauge_ids = graph[level]['deps']
        upgauge_int = [int(item) for item in upgauge_ids]
        # l.info(f"Updating the following upstream gauges: {upgauge_int}")
        
        # select matching row from random_params
        # masking matches exact param values, rather than indexing
        _mask = pd.Series([True] * len(random_params), index=random_params.index)
        
        for key, value in params.items():
            _mask = _mask & (random_params[key] == value)
        
        random_params_sel = random_params[_mask]
        
        # check if the matching is unique
        if len(random_params_sel) != 1:
            raise ValueError(f"Matching row in random_params is not unique! {len(random_params_sel)} rows found.")
        
        # for each upstream gauge, set random param from random_params_sel
        for upgauge in upgauge_int:
            vds_upgauge = vds[vds.value == upgauge]
            mask_up = ds.raster.geometry_mask(vds_upgauge)
            
            #ic(random_params_sel[str(float(upgauge))])
            #ic(random_params_sel[str(float(upgauge))].apply(eval))
            #ic(random_params_sel[str(float(upgauge))].iloc[0])
            #ic(ast.literal_eval(random_params_sel[str(float(upgauge))].iloc[0]))
            # select random paramset
            random_paramset = ast.literal_eval(random_params_sel[str(float(upgauge))].iloc[0])
            
            # modify staticmap based on the random paramset
            for key, value in random_paramset.items():
                da = ds[params_sname_to_lname[key]]
                
                if params_sname_to_method[key] == "mult":
                    da.values[mask_up] *= value
                    
                elif params_sname_to_method[key] == "set":
                    da.values[mask_up] = value
                
                elif params_sname_to_method[key] == "add":
                    da.values[mask_up] += value          
                
                ds[params_sname_to_lname[key]] = da
    
    ds.to_netcdf(out)
    l.info(f"Writing dataset to {out}")

    #make sure the lakes h-q rating relation is also copied
    for lin, lout in zip(lakes_in, lakes_out):
        if not os.path.exists(lout):
            shutil.copy(lin, lout)
            # l.info(f"Copying {lin} to {lout}")


if __name__ == "__main__":

    l=setup_logging('data/0-log', '03-set_calib_params.log')
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            # l.info(f"mod.params.params: {mod.params.params}")
            main(
                l,
                random_params=mod.input.random_params,
                p=mod.params.dataset,
                params=mod.params.params,
                params_sname=mod.params.params_sname,
                params_lname=mod.params.params_lname,
                params_method=mod.params.params_method,
                level=mod.params.level,
                graph=mod.params.graph,
                sub_catch=mod.params.sub_catch,
                lakes_in=mod.params.lake_in,
                lakes_out=mod.output.lake_out, 
                out=mod.output.staticmaps,
            )

        else:
            p = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticmaps/staticmaps.nc"
            p_out = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/staticmaps_test.nc"
            params_lname, params_method, all_level_df = create_set_all_levels(last_level=5, RECIPE="config/LHS_calib_recipe.json", N_SAMPLES=10, OPTIM='random-cd')
            sub_catch = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticgeoms/subcatch_Hall.geojson"
            lakes_in = ["/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticmaps/lake_hq_1.csv",
                        "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticmaps/lake_hq_2.csv",
                        "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticmaps/lake_hq_3.csv"]
            graph = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_levels_graph.json"))
            graph_pred = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_pred_graph.json"))
            graph_node = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_nodes_graph.json"))
            level = 'level1'
            best_params_previous = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/best_params.csv"
            params_df = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/paramspace.csv"
            random_params = Path('data/2-interim','calib_data', level, 'random_params.csv')
            params_df = pd.read_csv(params_df, index_col=0)
            # ic(params_df)
            # ic(params_df.iloc[0])
            params_cols = params_df.columns
            params_values = params_df.iloc[0].values
            params = dict(zip(params_cols, params_values))
            
            params_sname = list(all_level_df.columns)
            #replace 'nl' with 'nl,nf'
            params_sname = [param.replace('nl', 'nl,nf') for param in params_sname]
            params_sname = tuple(params_sname)
            lakes_out = [f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_1.csv",
                        f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_2.csv",
                        f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_3.csv"]
            
            main(
                l,
                p=p,
                params=params,
                params_sname=params_sname,
                params_lname=params_lname,
                params_method=params_method,
                random_params=random_params,
                level=level,
                graph=graph,
                sub_catch=sub_catch,
                lakes_in=lakes_in,
                lakes_out=lakes_out, 
                out=p_out,
            )
    except Exception as e:
        l.error(traceback.format_exc())
        raise e


