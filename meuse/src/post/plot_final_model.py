#TODO: organize functions to plot hydrograph, peak errors, etc.

from pathlib import Path

import numpy as np
import xarray as xr

from func_plot_signature import plot_signatures, plot_hydro
from typing import List

def main(
    md_data: Path | str,
    obs_data: Path | str,
    gauges: List[str],
    starttime: str,
    endtime: str,
    period_start: tuple | list,
    period_length: tuple | list,
    period_unit: str,
    output_dir: Path | str,
):


    """_summary_"""
    # cfg path
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()

    md_ds = xr.open_dataset(md_data)
    obs_ds = xr.open_dataset(obs_data)
    # obs_ds.Q.values = obs_ds.Q.values * (0.3048**3) #TODO: no need for this?

    # Convert gauges to integers for obs_ds
    gauges_int = [int(gauge) for gauge in gauges]

    temp = obs_ds.sel(wflow_id=gauges_int, time=slice(starttime, endtime), runs='Obs.')
    da = xr.DataArray(
        np.expand_dims(temp.Q.values, 2),
        coords = {
            "time": temp.time.values,
            "wflow_id": gauges_int,
            "runs": ["Obs.",]
        },
        dims = ["time", "wflow_id", "runs"]
    )
    temp = md_ds.sel(Q_gauges_Hall=gauges, time=slice(starttime, endtime))
    da = xr.concat(
        [
            da, 
            xr.DataArray(
                np.expand_dims(temp.Q.values, 2),
                coords = {
                    "time": temp.time.values,
                    "wflow_id": gauges_int,
                    "runs": ["Md.",]
                },
                dims = ["time", "wflow_id", "runs"],
            ),
        ],
        dim="runs",
    )
    ds = da.to_dataset(name="Q")
    
    return ds

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
    
    # #TODO: add a function to plot interactive hydrograph
    
    
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
    gauges = ['4', '801', '10', '11', '12', '16']
    starttime = '2008-08-01T00:00:00'
    endtime = '2018-02-22T00:00:00'
    period_start = ['2008-02-01', '2010-02-01', '2012-02-01']
    period_length = [365, 365, 365]
    period_unit = 'h'
    output_dir = work_dir / 'figures'
    
    ds = main(
        md_data=md_data,
        obs_data=obs_data,
        gauges=gauges,
        starttime=starttime,
        endtime=endtime,
        period_start=period_start,
        period_length=period_length,
        period_unit=period_unit,
        output_dir=output_dir,
    )
    
    from peak_timing import plot_peaks_ts
    import pandas as pd
    from datetime import datetime
    from metrics.run_peak_metrics import store_peak_info
    
    df_GaugeToPlot = pd.read_csv(work_dir / 'wflow_id_add_HBV_new.csv')
    start = datetime.strptime('2008-01-01', '%Y-%m-%d')
    end = datetime.strptime('2018-02-22', '%Y-%m-%d')
    
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
    run_keys = ds.runs.values
    color_dict = {f'{key}': color_list[i] for i, key in enumerate(run_keys)}
    
    
    peak_dict = store_peak_info(ds.sel(time=slice(start,end)), 'wflow_id', 72)
    
    
    plot_peaks_ts(ds,
                  df_GaugeToPlot,
                  start, 
                  end,
                  Folder_plots=output_dir,
                  color_dict=color_dict,
                  run_keys=run_keys,
                  peak_dict=peak_dict,
                  savefig=True,
                  id_key='wflow_id')
    

    
    
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
