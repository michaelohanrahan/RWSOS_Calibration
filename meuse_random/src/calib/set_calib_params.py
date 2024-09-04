
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

def handle_cs(params_sname, params_lname, params_method):
    '''
    Handle comma-separated items in params_sname and params_lname,
    adjusting params_method accordingly.
    '''
    new_sname = []
    new_lname = []
    new_method = []

    for sname, lname, method in zip(params_sname, params_lname, params_method):
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
    
    for key, value in params.items():
        # vds[key] = value   # what is this line used for?
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
        _mask = pd.Series([True] * len(random_params))
        
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
            
            # select random paramset
            random_paramset = random_params_sel[str(float(upgauge))].apply(eval)[0]
            
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
            l.info(f"Copying {lin} to {lout}")


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
            from snakemake.utils import Paramspace
            import json
            import pickle as pk
            
            work_dir = Path(r'p:\11209265-grade2023\wflow\UNREAL_TEST_DATA')
            csv = pd.read_csv(work_dir / 'best_10params.csv', index_col='gauges')
            csv.index=csv.index.astype(int)
            
            p = work_dir / 'staticmaps/staticmaps.nc'
            
            with open(work_dir/'create_set_params.pkl', 'rb') as f:
                d = pk.load(f)
            params_lname = d['lnames']
            params_method = d['methods']
            df = d['ds']
            
            paramspace = Paramspace(df)
            params = paramspace.instance
            
            level = 'level0'
            graph = json.load(open(Path(work_dir / 'Hall_levels_graph.json')))
            nodes = json.load(open(Path(work_dir / 'Hall_nodes_graph.json')))
            sub_catch = work_dir / 'subcatch_Hall.geojson'
            
            gauges = graph[level]['elements']
            pprint(gauges)
            csv['level'] = level
            # print(csv.head(1))
            def create_param(df):
                param = {col:round(random.uniform(0, 2), 2) for col in df.columns}
                return param
            
            data = {
                'gauges': list(gauges),
                **{'Top_{}'.format(i): str(create_param(df)) for i in range(10+1)},
            }
            
            best_params = pd.DataFrame(data)

            main(
                l,
                p=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\3-input\staticmaps\staticmaps.nc",
                params=params,
                params_lname=params_lname,
                params_method=params_method,
                random_params=random_params,
                level=mod.params.level,
                graph=mod.params.graph,
                sub_catch=mod.params.sub_catch,
                lakes_in=mod.params.lake_in,
                lakes_out=mod.output.lake_out, 
                out=mod.output.staticmaps,
            )
    except Exception as e:
        l.error(traceback.format_exc())
        raise e


