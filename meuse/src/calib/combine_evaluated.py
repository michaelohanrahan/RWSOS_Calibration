import xarray as xr
from pathlib import Path
from setuplog import setup_logging
import traceback
import pandas as pd
import timeit
from dask.diagnostics import ProgressBar
from tqdm import tqdm
from tqdm import trange
from tqdm.contrib.logging import logging_redirect_tqdm

def main(
    l,
    files: list,
    out_nc: str | Path,
    best_n: int = 10,
    method: str = 'ordinal',
):
    l.info(f"Combining {len(files)} files to {out_nc}")
    t0 = timeit.default_timer()
    ds = xr.open_mfdataset(files, combine="by_coords")
    l.info(f"mfdataset opened")
    l.info(f"Opening files took {timeit.default_timer() - t0:.2f} seconds")
    t1 = timeit.default_timer()
    l.info(f"Computing")
    with ProgressBar():
        ds = ds.compute()
    l.info(f"Computing took {timeit.default_timer() - t1:.2f} seconds")
    l.info(f"writing to {out_nc}")
    t2 = timeit.default_timer()
    with ProgressBar():
        ds.to_zarr(Path(out_nc), mode='w')
    l.info(f"Writing to zarr took {timeit.default_timer() - t2:.2f} seconds")
    l.info(f"Combining files took {timeit.default_timer() - t0:.2f} seconds")

    # retries = 5 
    # for attempt in range(retries):
    #     try:
            # with xr.open_dataset(out_nc) as ds:
                # Automatically detect parameter dimensions by excluding known non-parameter dimensions
    non_param_dims = ['gauges']
    param_dims = [dim for dim in ds.dims if dim not in non_param_dims]
    
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
    with open(Path(out_nc).parent / "eval.done", "w") as f:
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
            main(l,
                ['p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.0/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.0/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.0/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.0/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.0/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.0/rd~1.75/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~0.25/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.0/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.0/f~1.5/rd~1.75/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~0.25/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.0/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.2/rd~1.75/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~0.25/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.0/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.5/rd~1.75/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~0.25/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.0/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~0.75/rd~1.75/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~0.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.0/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~0.25/st~1.5/nr~1.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~0.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~0.5/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~1.0/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~1.0/ml~0.3/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~1.5/ml~0.0/performance.nc',
 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~1.5/f~1.0/rd~1.0/st~0.5/nr~1.5/ml~0.3/performance.nc'],
                out_nc=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\performance_test.zarr",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise