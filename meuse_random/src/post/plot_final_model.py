#TODO: add benchmarks to main function

from pathlib import Path
from typing import List
from icecream import ic
from glob import glob

import numpy as np
import xarray as xr
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
# self defined functions
from peak_timing import plot_peaks_ts, peak_timing_for_runs, plot_peak_timing_distribution
from metrics.run_peak_metrics import store_peak_info
from func_plot_signature import plot_signatures


def main(
    md_data: List[Path | str],
    obs_data: Path | str,
    starttime: str,
    endtime: str,
    GaugeToPlot: Path | str, # path to the dataframe of gauges to plot (wflow_id_add_HBV_new.csv)
    savefig: bool,
    plotfig: bool,
    output_dir: Path | str,
):

    """
    Main function to process and plot hydrological data.
    Args:
        md_data: List[xr.Dataset] A list of nc paths containing results from model runs.
        obs_data (Path | str): Path to the observed data file.
        starttime (str): Start time for data processing and plotting.
        endtime (str): End time for data processing and plotting.
        GaugeToPlot (Path | str): Path to the CSV file containing gauge information.
        savefig (bool): Flag to save generated figures.
        plotfig (bool): Flag to display generated figures.
        output_dir (Path | str): Directory to save output files and figures.
    Returns:
        None.

    This function reads model and observed data, processes it for specified gauges and time period,
    and prepares a combined dataset. It can optionally generate and save various plots including
    hydrographs, peak timing analysis, and hydrological signatures.
    """
    # cfg path
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()
    
    
    md_ds = []
    level_list = []
    
    for file in md_data:
        try:
            md_ds.append(xr.open_dataset(file))
            level_list.append(Path(Path(file).parent).parent.name.split('_')[-2])
        except:
            print(f'{file} is not a valid file')
            continue
    
    level_list = [level.replace('level-1', 'base') for level in level_list]
    
    obs_ds = xr.open_dataset(obs_data)
    
    df_GaugeToPlot = pd.read_csv(GaugeToPlot)
    # gauges = list(df_GaugeToPlot['wflow_id'].values)  # integers
    gauges=[16]
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
            ) for file, level in zip(md_ds, level_list)
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
    
    
    # plot interactive hydrograph
    color_list = ['#377eb8', 
              '#ff7f00', 
              '#4daf4a', 
              '#f781bf', 
              '#a65628', 
              '#984ea3', 
              '#999999', 
              '#e41a1c', 
              '#dede00', 
              '#ff7f00', 
              '#a65628', 
              '#f781bf'
              ]
    
    start = datetime.strptime(starttime, '%Y-%m-%d')
    end = datetime.strptime(endtime, '%Y-%m-%d')
    run_keys = ds.runs.values
    color_dict = {f'{key}': color_list[i] for i, key in enumerate(run_keys)}
    peak_dict = store_peak_info(ds.sel(time=slice(starttime,endtime)), 'wflow_id', 72)
    

    plot_peaks_ts(ds,
                  df_GaugeToPlot,
                  start, 
                  end,
                  Folder_plots=output_dir,
                  color_dict=color_dict,
                  run_keys=run_keys,
                  peak_dict=peak_dict,
                  savefig=savefig,
                  id_key='wflow_id')
    
    # plot peak timing error
    peak_timing_for_runs(ds,
                         df_GaugeToPlot, 
                         output_dir, 
                         peak_dict=peak_dict,
                         plotfig=plotfig,
                         savefig=savefig)
    
    plot_peak_timing_distribution(ds,
                                  df_GaugeToPlot,
                                  run_keys,
                                  peak_dict,
                                  color_dict,
                                  output_dir,
                                  cumulative=False,
                                  savefig=True,
                                  id_key='wflow_id')
    plot_peak_timing_distribution(ds,
                                  df_GaugeToPlot,
                                  run_keys,
                                  peak_dict,
                                  color_dict,
                                  output_dir,
                                  cumulative=True,
                                  savefig=True,
                                  id_key='wflow_id')

if __name__ == "__main__":
    
    work_dir = Path(r'p:\11209265-grade2023\wflow\wflow_meuse_julia').as_posix()
    files = glob(str(Path(work_dir, 'best_run_level*_result', 'output_run', 'output_scalar.nc')))
    # files = [file for file in files if 'level0' not in file]
    
    obs_data = Path(r"P:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random\data\1-external\discharge_hourlyobs_smoothed.nc")
    GaugeToPlot = Path(work_dir, 'best_run_level-1_result','wflow_id_add_HBV_new.csv')
    starttime = '2008-08-01'
    endtime = '2018-02-22'
    output_dir = Path(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\5-visualization\best_params_smoothed").as_posix()
    
    ds = main(
        md_data=files,
        obs_data=obs_data,
        starttime=starttime,
        endtime=endtime,
        GaugeToPlot=GaugeToPlot,
        savefig=True,
        plotfig=True,
        output_dir=output_dir,
    )
    