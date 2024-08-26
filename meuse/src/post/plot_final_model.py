#TODO: add benchmarks to main function

from pathlib import Path
from typing import List

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
    md_data: Path | str,
    obs_data: Path | str,
    starttime: str,
    endtime: str,
    GaugeToPlot: Path | str, # path to the dataframe of gauges to plot (wflow_id_add_HBV_new.csv)
    savefig: bool,
    plotfig: bool,
    # period_start: tuple | list,
    # period_length: tuple | list,
    # period_unit: str,
    output_dir: Path | str,
):

    """
    Main function to process and plot hydrological data.
    Args:
        md_data (Path | str): Path to the model output data file.
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

    md_ds = xr.open_dataset(md_data)
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
    
    # prepare md data-array
    temp = md_ds.sel(Q_gauges_Hall=gauges_str, time=slice(starttime, endtime))
    da = xr.concat(
        [
            da, 
            xr.DataArray(
                np.expand_dims(temp.Q.values, 2),
                coords = {
                    "time": temp.time.values,
                    "wflow_id": gauges,
                    "runs": ["Md.",]
                },
                dims = ["time", "wflow_id", "runs"],
            ),
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
              '#f781bf']
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
    
    # not fully working
    plot_peak_timing_distribution(run_keys, peak_dict, color_dict, output_dir)

    # hydro_periods = []
    # for idx, _s in enumerate(period_start):
    #     _s = np.datetime64(_s)
    #     _e = _s + np.timedelta64(period_length[idx], period_unit)
    #     hydro_periods.append(
    #         (
    #             np.datetime_as_string(_s, unit="h"),
    #             np.datetime_as_string(_e, unit="h"),
    #         )
    #     )
    
    # pass
    # ### prepare dataset to make plots
    
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

    # # ##%%
    # plot_colors = colors[:1]
    # for station_id in gauges:
    #     # Extract station name
    #     try:
    #         station_name = str(station_id)
    #     except KeyError:
    #         station_name = str(station_id)

    #     # Get data for that stations
    #     dsq = ds.sel(stations=station_id)
    #     # Prep labels
    #     labels = np.delete(dsq.runs.values, np.where(dsq.runs.values == "Obs."))

    #     # Plot hydrograpgh
    #     plot_hydro(
    #         dsq=dsq,
    #         start_long=starttime.split("T")[0],
    #         end_long=endtime.split("T")[0],
    #         start_1=hydro_periods[0][0],
    #         end_1=hydro_periods[0][1],
    #         start_2=hydro_periods[1][0],
    #         end_2=hydro_periods[1][1],
    #         start_3=hydro_periods[2][0],
    #         end_3=hydro_periods[2][1],
    #         labels=labels,
    #         colors=plot_colors,
    #         Folder_out=output_dir,
    #         station_name=station_name,
    #         station_id=station_id,
    #         save=True,
    #     )

    #     # plot signatures
    #     plot_signatures(
    #         dsq=dsq,
    #         labels=labels,
    #         colors=plot_colors,
    #         Folder_out=output_dir,
    #         station_name=station_name,
    #         station_id=station_id,
    #         save=True,
    #     )

if __name__ == "__main__":
    
    work_dir = Path(r'p:\11209265-grade2023\wflow\wflow_meuse_julia\best_run_level0_result')
    md_data = work_dir / 'output_run/output_scalar.nc'
    obs_data = work_dir / 'discharge_hourlyobs_HBV_combined.nc'
    GaugeToPlot = work_dir / 'wflow_id_add_HBV_new.csv'
    starttime = '2008-08-01'
    endtime = '2018-02-22'
    output_dir = work_dir / 'figures'
    
    ds = main(
        md_data=md_data,
        obs_data=obs_data,
        starttime=starttime,
        endtime=endtime,
        GaugeToPlot=GaugeToPlot,
        savefig=False,
        plotfig=False,
        output_dir=output_dir,
    )
    
    
    # if "snakemake" in globals():
    #     mod = globals()["snakemake"]
        
    #     main(
    #         mod.input.scalar,
    #         mod.params.observed_data,
    #         mod.params.gauges,
    #         mod.params.starttime,
    #         mod.params.endtime,
    #         mod.params.period_startdate,
    #         mod.params.period_length,
    #         mod.params.period_unit,
    #         mod.params.output_dir,
    #     )
    
    # else:
    #     main(
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/Intermittent_result/run_default/output_scalar.nc",
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/data/GAUGES/discharge_obs_combined.nc",
    #         ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000', '12106700', '12115000', '12121600'],
    #         "2012-01-01T00:00:00",
    #         "2018-12-31T00:00:00",
    #         ["2012-09-01", "2013-09-01", "2015-09-01"],
    #         [365, 365, 365],
    #         "D",
    #         "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/Intermittent_result/figures"
    #     )
