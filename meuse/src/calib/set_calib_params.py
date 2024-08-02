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

def create_coords_from_params(params: dict):
    """Create coordinate MultiIndex from parameter dictionary."""
    keys = list(params[0].keys())
    params_values = list(zip(*[_dict.values() for _dict in params]))
    coord = pd.MultiIndex.from_arrays(params_values, names=keys)
    return coord

def fake_performance_nc():
    params = [
        {'param1': 1, 'param2': 21, 'param3': 0.1},
        {'param1': 2, 'param2': 6, 'param3': 0.2},
        {'param1': 3, 'param2': 12, 'param3': 0.3},
        {'param1': 4, 'param2': 18, 'param3': 0.4},
        {'param1': 5, 'param2': 24, 'param3': 0.5},
        {'param1': 6, 'param2': 30, 'param3': 0.6},
        {'param1': 7, 'param2': 36, 'param3': 0.7},
        {'param1': 8, 'param2': 42, 'param3': 0.8},
        {'param1': 9, 'param2': 48, 'param3': 0.9},
        {'param1': 10, 'param2': 54, 'param3': 1.0}
    ]

    param_coords = create_coords_from_params(params)
    gauges = ['gauge1', 'gauge2']
    metrics = ['kge', 'rld']

    # Initialize metric_values with the correct shape
    metric_values = {
        'kge': np.random.rand(len(params), len(gauges)),
        'rld': np.random.rand(len(params), len(gauges))
    }

    res = np.random.rand(len(params), len(gauges))
    weights = [0.5, 0.5]

    ds = None
    for metric in metrics:
        da = xr.DataArray(
            metric_values[metric], 
            coords=[('params', param_coords), ('gauges', gauges)], 
            dims=['params', 'gauges']
        )
        da.name = metric
        if ds is None:
            ds = da.to_dataset()
            continue

        ds[metric] = da

    ds["euclidean"] = xr.DataArray(
        res, 
        coords=[('params', param_coords), ('gauges', gauges)], 
        dims=['params', 'gauges']
    )

    ds = ds.assign_attrs(
        {
            "metrics": metrics,
            "weights": weights,
        }
    )

    # Fill NaNs with a default value before unstacking
    ds = ds.fillna(0)
    ds = ds.unstack()

    return ds

# Example usage
ds = fake_performance_nc()
print(ds)

def upper_level_random_params(l, 
                              ds, 
                              best_params, 
                              level, 
                              graph, 
                              vds, 
                              params_lnames,
                              params_method,
                              test=False)->xr.Dataset:
    '''
    if upper level is present, randomly select some paramters from the best_params
    we will have to save to the outfolder a csv of params per gauge_id chosen
    '''
    if not os.path.exists(best_params) or level == 'level0':
        l.info(f"Best params file not found or level is level0, skipping upper level random params")
        return ds
    else:
        # if not os.path.exists(best_params):
        #     l.error(f"Best params file not found: {best_params}")
        #     raise FileNotFoundError(f"Best params file not found: {best_params}")
        best_params = fake_performance_nc()
        l.info(f"Loading best params from {best_params}")
        
        #load the best params ()
        best_params_ds = xr.open_dataset(best_params)
        
        #approximated data input
        # <xarray.Dataset> Size: 48kB
        # Dimensions:    (param1: 10, param2: 10, param3: 10, gauges: 2)
        # Coordinates:
        # * param1     (param1) int64 80B 1 2 3 4 5 6 7 8 9 10
        # * param2     (param2) int64 80B 6 12 18 21 24 30 36 42 48 54
        # * param3     (param3) float64 80B 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
        # * gauges     (gauges) <U6 48B 'gauge1' 'gauge2'
        # Data variables:
        #     kge        (gauges, param1, param2, param3) float64 16kB nan nan ... 0.8428
        #     rld        (gauges, param1, param2, param3) float64 16kB nan nan ... 0.3971
        #     euclidean  (gauges, param1, param2, param3) float64 16kB nan nan ... 0.1172
        # Attributes:
        #     metrics:  ['kge', 'rld']
        #     weights:  [0.5, 0.5]
        
        # rename level to level{-1}
        lsplit = int(level.split('level')[0])
        level = f'level{lsplit-1}'
        
        
        gauge_id = graph[level]["elements"]
        gauge_int = [
            int(item) for item in gauge_id
        ]
        l.info(f"Updating the following upstream gauges: {gauge_int}")
        vds_u = vds[vds.value.isin(gauge_int)]
        
        #random int
        random_ints = [random.randint(1, 10) for _ in gauge_int]
        gauge_dict = {gauge: num for gauge, num in zip(gauge_int, random_ints)}
        col_method = {col: method for col, method in zip(params_lnames, params_method)}
        
        #TODO: finish
        for gauge, index in gauge_dict.items():
            #only target one gauge_id at a time
            mask = ds.raster.geometry_mask(vds_u[vds_u.value == gauge])
            param_da = best_params_ds["euclidean"].sel(gauges=gauge)
            #what are the associated params for each euc value
            
            df = #indexed by euc value and params name in cols
            #sort and retain the highest euc values
            df_sorted = df.sort_values(by=gauge, ascending=False)
            params = df_sorted.iloc[index]
            params = params.to_dict()
            
            for idx, (key, value) in enumerate(params.items()):
                da = ds[params_lnames[idx]]
                if col_method[idx] == "mult":
                    da.values[mask] *= value
                elif col_method[idx] == "set":
                    da.values[mask] = value
                elif col_method[idx] == "add":
                    da.values[mask] += value          
                ds[params_lnames[idx]] = da
            


def main(
    l,
    best_params: Path,
    p: Path | str,
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    level: str,
    graph: dict,
    # gauge_ids: tuple | list,
    sub_catch: Path | str,
    lakes_in: Path | str,
    lakes_out: Path | str,
    out: Path | str,
):
    """
    Load the dataset (staticmaps) then update params for the desired gauge_ids
    the paramslname should be the same as the dataset key
    """
    gauge_ids_d = graph[level]["elements"]
    
    # Load the current main dataset
    ds = xr.open_dataset(p)
    l.info(f"Loading dataset from {p}")
    
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
    
    # Create the data mask
    mask = ds.raster.geometry_mask(vds_d)

    # Loop through the parameters
    for idx, (key, value) in enumerate(params.items()):
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

    ds = upper_level_random_params(l,
                                   ds, 
                                   best_params, 
                                   level, 
                                   graph, 
                                   vds,
                                   params_lname,
                                   params_method,
                                   test=False)
    
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

