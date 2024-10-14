#TODO: add benchmarks to main function
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
    run_list: List[str],
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
    
    for run, file in zip(run_list, md_data):
        try:
            md_ds.append(xr.open_dataset(file))
            print(f'"{run}" run added via file:\n{file}')
        except:
            print(f'{file} is not a valid file')
            continue
    
    obs_ds = xr.open_dataset(obs_data)
    
    df_GaugeToPlot = pd.read_csv(GaugeToPlot)
    gauges = list(df_GaugeToPlot['wflow_id'].values)  # integers
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
              '#f781bf',
              '#0072b2',
              ]
    
    start = datetime.strptime(starttime, '%Y-%m-%d')
    end = datetime.strptime(endtime, '%Y-%m-%d')
    run_keys = ds.runs.values
    color_dict = {f'{key}': color_list[i] for i, key in enumerate(run_keys)}
    peak_dict = store_peak_info(ds.sel(time=slice(starttime,endtime)), 'wflow_id', 72)

    # plot_peaks_ts(ds,
    #               df_GaugeToPlot,
    #               start, 
    #               end,
    #               Folder_plots=output_dir,
    #               color_dict=color_dict,
    #               run_keys=run_keys,
    #               peak_dict=peak_dict,
    #               savefig=savefig,
    #               id_key='wflow_id')
    
    # # plot peak timing error
    # peak_timing_for_runs(ds,
    #                      df_GaugeToPlot, 
    #                      output_dir, 
    #                      peak_dict=peak_dict,
    #                      plotfig=plotfig,
    #                      savefig=savefig)
    
    plot_peak_timing_distribution(ds,
                                  df_GaugeToPlot,
                                  run_keys,
                                  peak_dict,
                                  color_dict,
                                  output_dir,
                                  cumulative=True,
                                  savefig=True,
                                  id_key='wflow_id')
    
    # # plot signature
    # colors = [
    #     # "#a6cee3",
    #     # "#1f78b4",
    #     # "#b2df8a",
    #     # "#33a02c",
    #     # "orange",
    #     # "#fb9a99",
    #     # "#e31a1c",
    #     # "#fdbf6f",
    #     "#ff7f00",
    #     "#cab2d6",
    #     "#6a3d9a",
    #     "#ffff99",
    #     "#b15928",
    # ]
    # plot_colors = colors[:1]
    
    # for station_id in gauges:
    #     # Extract station name
    #     station_name = df_GaugeToPlot.loc[df_GaugeToPlot['wflow_id']==station_id, 'location'].values[0]

    #     # Get data for that stations
    #     dsq = ds.sel(wflow_id=station_id)
    #     # Prep labels
    #     labels = np.delete(dsq.runs.values, np.where(dsq.runs.values == "Obs."))

    #     plot_signatures(
    #         dsq=dsq,
    #         labels=labels,
    #         colors=plot_colors,
    #         Folder_out=output_dir,
    #         station_name=station_name,
    #         station_id=station_id,
    #         save=True,
    #     )
    
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
    
    ds = main(
        md_data=files,
        obs_data=obs_data,
        starttime=starttime,
        endtime=endtime,
        run_list=run_list,
        GaugeToPlot=GaugeToPlot,
        savefig=True,
        plotfig=True,
        output_dir=output_dir,
    )
    
    
    
    ds.to_netcdf(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider\data\4-output\combined_output.nc')