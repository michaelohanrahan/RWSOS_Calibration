import xarray as xr
from pathlib import Path
from setuplog import setup_logging
import traceback
import pandas as pd
import dask
from dask.diagnostics import ProgressBar
from glob import glob
from icecream import ic
import os 
import timeit
import numpy as np
import random

def calculate_euclidean_threshold(l, file_list, percentile=10):
    """
    Opens a multi-file dataset, selects a random subset of files,
    and calculates the 10th percentile of the 'euclidean' values.

    Parameters:
    - file_list: list of str
        List of file paths to open as a dataset.
    - percentile: int
        The percentile to calculate (default is 10).

    Returns:
    - threshold: float
        The nth percentile threshold value of the 'euclidean' values.
    """
    # Randomly select a subset of files
    selected_files = random.sample(file_list, len(file_list))  # Use all files in the list

    l.info(ic(len(selected_files)))
    # Open the selected files as a single xarray dataset
    ds = xr.open_mfdataset(selected_files, combine='by_coords')
    l.info(ic(ds))
    # Extract the 'euclidean' values
    euclidean_values = ds['euclidean'].values.flatten()  # Flatten to 1D array

    # Calculate the specified percentile
    threshold = np.nanpercentile(euclidean_values, percentile)

    return threshold

def del_existing(out_nc):
    # if not os.path.exists(out_nc):
    if os.path.exists(out_nc):
       if os.path.isfile(out_nc):
           os.remove(out_nc)  # Remove the existing file
           l.info(f"Removed existing file: {out_nc}")
       elif os.path.isdir(out_nc):
           l.info(f"Output path already exists as a directory: {out_nc}")
           
def select_random_files(file_list, num_files=1000):
    """
    Selects a random number of files from the provided list.

    Parameters:
    - file_list: list of str
        List of file paths to select from.
    - num_files: int
        The number of random files to select.

    Returns:
    - selected_files: list of str
        A list of randomly selected file paths.
    """
    # Ensure num_files does not exceed the length of file_list
    num_files = min(num_files, len(file_list))
    
    # Randomly select the specified number of files
    selected_files = random.sample(file_list, num_files)
    
    return selected_files

def main(
    l,
    files: list,
    out_nc: str|Path,
    best_n: int = 10,
    n_subsample:int=100,
    sub_pct:int=10
):
    for file in files:
        del_existing(file)
    # random_files = select_random_files(files, num_files=n_subsample)
    
    # s_euc = timeit.default_timer()
    # euc_threshold = calculate_euclidean_threshold(l, random_files, percentile=sub_pct)
    # e_euc = timeit.default_timer()
    # l.info(f"calculated euclidean threshold to be: {euc_threshold}, from a subsample of {n_subsample} indicating a percentile of {sub_pct}")
    # l.info(f"calculation time: {e_euc-s_euc}s")
    

    l.info(f"Combining {len(files)} files to {out_nc}")
    t0 = timeit.default_timer()
    chunksize = {'gauges': -1, 'ksat': 1, 'f': 1, 'rd': 1, 'st': 1, 'nr': 1, 'ml': 1, 'nl': 1, 'nf': 1}
    ds = xr.open_mfdataset(files, combine="by_coords", chunks=chunksize)
    l.info(f"mfdataset opened")
    l.info(ic(ds))
    l.info(ic(ds['euclidean']))
    a
    l.info(f"Opening files took {timeit.default_timer() - t0:.2f} seconds")
    t1 = timeit.default_timer()
    l.info(f"Computing")
    
    # with ProgressBar():
    #     ds = ds.compute()
    
    l.info(f"Computing took {timeit.default_timer() - t1:.2f} seconds")
    l.info(f"writing to {out_nc}")
    t2 = timeit.default_timer()
    
    with ProgressBar():
        ds.to_zarr(Path(out_nc), mode='w')
    
    l.info(f"Writing to zarr took {timeit.default_timer() - t2:.2f} seconds")
    l.info(f"Combining files took {timeit.default_timer() - t0:.2f} seconds") 
    
    # if os.path.exists(out_nc):
    #     l.info(f"opened ")
    #     ds = xr.open_zarr(out_nc)
    
    # Automatically detect parameter dimensions by excluding known non-parameter dimensions
    non_param_dims = ['gauges']
    param_dims = [dim for dim in ds.dims if dim not in non_param_dims]
    # print(param_dims)
    
    # Flatten the parameter dimensions into a single dimension
    flattened = ds.stack(params=param_dims).compute()
    print(flattened)
    # Rank the 'euclidean' values per gauge
    ranked = flattened['euclidean'].rank(dim='params')

    best_params = []
    for gauge in ds['gauges'].values:
        # Drop NaN values for this specific gauge to retain only valid parameter sets
        valid_euclidean = flattened['euclidean'].sel(gauges=gauge).dropna(dim='params')
        valid_ranked = ranked.sel(gauges=gauge).dropna(dim='params')
        top_n = min(best_n, len(valid_euclidean))
        
        best_for_gauge = {'gauge': gauge}
        for idx in range(top_n):
            param_vals = valid_ranked.where(valid_ranked == idx + 1, drop=True).params
            param_dict = {param: param_vals[param].values[0] for param in param_dims}
            best_for_gauge[f'Top_{idx+1}'] = param_dict
        best_params.append(best_for_gauge)
    
    # Flatten the parameter dimensions into a single dimension
    flattened = ds.stack(params=param_dims)
    best_params = []
    with logging_redirect_tqdm(loggers=[l]):
        for gauge in tqdm(ds['gauges'].values, 
                        desc='Finding best parameters for each gauge',
                        total=len(ds['gauges'].values),
                        unit='gauge',
                        ):
            # Drop NaN values and sort by 'euclidean' for the current gauge
            sort = flattened['euclidean'].sel(gauges=gauge).dropna(dim='params').sortby(flattened['euclidean'].sel(gauges=gauge))
            #add preferred ranking method
            
            top_n = min(best_n, len(sort.values))
            best_for_gauge = {'gauge': gauge}
            for idx in range(top_n):
                param_vals = sort[idx + 1].params
                if param_vals.size > 0:
                    param_dict = {param: float(param_vals[param].values) for param in param_dims}
                    best_for_gauge[f'Top_{idx + 1}'] = param_dict
            
            best_params.append(best_for_gauge)
    #TODO: Add the level column to best_params csv
    # Convert the results to a DataFrame and save as CSV
    final_df = pd.DataFrame(best_params)
    folder = Path(out_nc).parent
    out = Path(folder, "best_params").with_suffix('.csv')
    final_df.to_csv(out, index=False)
    
    # Create an 'eval.done' file to indicate completion
    with open(Path(out_nc).parent / "level.done", "w") as f:
        f.write("")
    l.info(f"Saved {out_nc}")


if __name__ == "__main__":
    l = setup_logging("data/0-log", "05-combine_evaluated.log")
    try:
        if 'snakemake' in globals():
            main(
                l,
                snakemake.input.performance_files,
                snakemake.output.performance,
            )
        else:
            base_dir = r"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0"
            l.info(f"Base dir exists: {Path(base_dir).exists()}\n{base_dir}")
            
            specific_file = 'performance.nc'
            found = []
            for root, dirs, files in os.walk(base_dir):
                if specific_file in files:
                    found.append(f"{root}/{specific_file}")
            l.info(f"{len(found)} files found")
            main(l,
                found[:10],
                out_nc="/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/performance_test.zarr",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise