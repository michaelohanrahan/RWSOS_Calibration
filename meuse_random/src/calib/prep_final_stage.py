#TODO: to be tested!

from pathlib import Path
import pandas as pd
import geopandas as gpd
import xarray as xr
import shutil
import os
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
    p: Path | str,
    best_params: pd.DataFrame,
    topx,
    params_sname: tuple | list,
    params_lname: tuple | list,
    params_method: tuple | list,
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
    params = best_params.loc[('level5'), topx]
    
    # Load the geometries
    vds = gpd.read_file(sub_catch)
    vds = vds.astype({"value": int})
    
    # Set param for all gauges
    gauge_ids = params.index.values  # float
    gauge_int = [int(item) for item in gauge_ids] # int
    
    for gauge in gauge_int:
        vds_gauge = vds[vds.value == gauge]
        mask_up = ds.raster.geometry_mask(vds_gauge)
        
        # select random paramset
        gauge_paramset = ast.literal_eval(params[float(gauge)])
        
        # modify staticmap based on the random paramset
        for key, value in gauge_paramset.items():
            da = ds[params_sname_to_lname[key]]
            
            if params_sname_to_method[key] == "mult":
                da.values[mask_up] *= value
                
            elif params_sname_to_method[key] == "set":
                da.values[mask_up] = value
            
            elif params_sname_to_method[key] == "add":
                da.values[mask_up] += value          
            
            ds[params_sname_to_lname[key]] = da
    
    ds.to_netcdf(out)
    
    #make sure the lakes h-q rating relation is also copied
    for lin, lout in zip(lakes_in, lakes_out):
        if not os.path.exists(lout):
            shutil.copy(lin, lout)



if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        main(
            p=mod.params.dataset,
            best_params=mod.params.best_params,
            topx = lambda wildcards: wildcards.Topx,
            params_sname=mod.params.params_sname,
            params_lname=mod.params.params_lname,
            params_method=mod.params.params_method,
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
       
        params_sname = list(all_level_df.columns)
        #replace 'nl' with 'nl,nf'
        params_sname = [param.replace('nl', 'nl,nf') for param in params_sname]
        params_sname = tuple(params_sname)
        lakes_out = [f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_1.csv",
                    f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_2.csv",
                    f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/ksat~{params['ksat']}/f~{params['f']}/rd~{params['rd']}/st~{params['st']}/nr~{params['nr']}/ml~{params['ml']}/nl~{params['nl']}/nf~{params['nf']}/lake_hq_3.csv"]
        
        best_params_fn = Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider\data\2-interim\calib_data\level5\best_params.csv')
        best_params = pd.read_csv(best_params_fn, index_col=['level', 'gauge'])
        
        main(
            p=p,
            best_params=best_params,
            topx = 'Top1',
            params_sname=params_sname,
            params_lname=params_lname,
            params_method=params_method,
            sub_catch=sub_catch,
            lakes_in=lakes_in,
            lakes_out=lakes_out, 
            out=p_out,
        )