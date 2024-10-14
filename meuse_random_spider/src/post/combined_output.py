import os
from pathlib import Path
from typing import List
# from icecream import ic
from glob import glob
import re
import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt

def combined_output(
    md_data: List[Path | str],
    obs_data: Path | str,
    starttime: str,
    endtime: str,
    GaugeToPlot: Path | str | None, # path to the dataframe of gauges to plot (wflow_id_add_HBV_new.csv)
    run_list: List[str],
    output_dir: Path | str,
):

    """
    This function reads model (Topx, base model, hbv) and observed data, processes it for specified gauges and time period,
    and prepares a combined dataset.
    Args:
        md_data: List[xr.Dataset] A list of nc paths containing results from model runs.
        obs_data (Path | str): Path to the observed data file.
        starttime (str): Start time for data processing and plotting.
        endtime (str): End time for data processing and plotting.
        GaugeToPlot (Path | str): Path to the CSV file containing gauge information.
        output_dir (Path | str): Directory to save output files and figures.
    Returns:
        xr.Dataset: A dataset containing the combined output data.
    """
    # cfg path
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    
    
    md_ds = []
    
    for run, file in zip(run_list, md_data):
        try:
            md_ds.append(xr.open_dataset(file))
            print(f'"{run}" run added via file:\n{file}')
        except:
            print(f'{file} is not a valid file')
            continue
    
    obs_ds = xr.open_dataset(obs_data)
    
    if GaugeToPlot is not None:
        df_GaugeToPlot = pd.read_csv(GaugeToPlot)
        gauges = list(df_GaugeToPlot['wflow_id'].values)  # integers
    else:
        gauges = [int(gauge) for gauge in md_ds[0].Q_gauges_Hall.values] # integers
        
    gauges_str = [str(gauge) for gauge in gauges]  # Convert gauges to strings for md_ds
        
    # prepare obs data-array
    temp = obs_ds.sel(wflow_id=gauges, time=slice(starttime, endtime), runs='Obs.')
    da = xr.DataArray(
        np.expand_dims(temp.Q.values, 2),
        coords = {
            "time": temp.time.values,
            "wflow_id": gauges,
            "runs": ["Obs.",]
        },
        dims = ["time", "wflow_id", "runs"]
    )
    
    # prepare HBV data-array
    temp = obs_ds.sel(wflow_id=gauges, time=slice(starttime, endtime), runs='HBV')
    
    da = xr.concat(
        [
            da, 
            xr.DataArray(
                np.expand_dims(temp.Q.values, 2),
        coords = {
            "time": temp.time.values,
            "wflow_id": gauges,
            "runs": ["HBV",]
                },
                dims = ["time", "wflow_id", "runs"]
            ),
        ],
        dim="runs"
    )
    
    # prepare md data-array
    md_data_dict = {level: xr.DataArray(
                np.expand_dims(file.sel(Q_gauges_Hall=gauges_str,
                                        time=slice(starttime, endtime)
                                        ).Q.values, 2),
                coords = {
                    "time": file.sel(Q_gauges_Hall=gauges_str,
                                        time=slice(starttime, endtime)
                                        ).time.values,
                    "wflow_id": gauges,
                    "runs": [level,]
                },
                        dims = [
                            "time", 
                            "wflow_id", 
                            "runs"]
            ) for file, level in zip(md_ds, run_list)
                    }
    da = xr.concat(
        [
            da, 
            *md_data_dict.values()
            
        ],
        dim="runs",
    )
    
    # convert data-array to dataset
    ds = da.to_dataset(name="Q")
    
    return ds


def natural_sort_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]


if __name__ == "__main__":
    
    work_dir = Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider').as_posix()
    
    files = sorted(glob(str(Path(work_dir, 'data', '4-output', 'output_*', 'output_scalar.nc'))),
                                 key=natural_sort_key)

    run_list = ['_'.join(Path(top).parts[-2].split('_')[-2:]) for top in files]
    
    obs_data = Path(work_dir, 'data/1-external/discharge_hourlyobs_smoothed.nc')
    GaugeToPlot = Path(work_dir,'data', '4-output','wflow_id_add_HBV_new.csv')

    starttime = '2005-08-01'
    endtime = '2018-02-22'
    output_dir = Path(work_dir, 'data/5-visualization/best_run').as_posix()
    os.makedirs(output_dir, exist_ok=True)
    
    ds = combined_output(
        md_data=files,
        obs_data=obs_data,
        starttime=starttime,
        endtime=endtime,
        run_list=run_list,
        GaugeToPlot=None,
        # GaugeToPlot=GaugeToPlot,
        output_dir=output_dir,
    )
    
    
    ds.to_netcdf(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider\data\4-output\combined_output.nc')