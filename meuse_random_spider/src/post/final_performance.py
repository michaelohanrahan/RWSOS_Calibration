from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from timeit import default_timer as timer
from metrics.metrics import kge, nse, nse_log, nselog_mm7q, mae_peak_timing, mape_peak_magnitude
from metrics.metrics import _obs_peaks, _sim_peaks

def final_performance(
    ds: Path | str,
    gid: str,
    starttime: str,
    endtime: str,
    gauges: tuple | list,
    metrics: tuple | list,
    dry_month: list,
    window: int,
    output_dir: Path | str,
    filename: str,
    ):
    
    # OPEN THE DATASET AND SELECT THE TIME PERIOD
    ds = xr.open_dataset(ds)
    ds = ds.sel(time=slice(starttime, endtime))
    
    # get metrics
    METRICS = {metric: globals()[metric] for metric in metrics}
    
    # make sure the gauges are int
    gauges = [int(g) for g in gauges]
    
    # CALCULATE PEAKS
    if any("peak" in metric for metric in metrics):
        start = timer()
        obs_peaks = {
            g: _obs_peaks(ds.sel(runs='Obs.',wflow_id=g).Q)
            for g in gauges
        }
        peaks = {
            run: {
                g: _sim_peaks(sim=ds.sel({gid: g}).sel(runs=run).Q, 
                              obs=obs_peaks[g],
                              window=window)
                for g in gauges
            }
            for run in ds.runs.values if run != 'Obs.'
        }
        end = timer()
    else:
        peaks = None
    print(f"Peaks calculated in {end - start} seconds")
    
    # CALCULATE METRIC VALUES FOR EACH RUN AND STORE RESULTS FOR ALL RUNS
    all_run_results = {}
    
    for run in ds.runs.values:
        if run == 'Obs.':
            continue  # Skip the 'Obs.' run as it's the observation data
        md = ds.sel(runs=run)
        obs = ds.sel(runs='Obs.')
        
        run_results = {}
        
        for metric in metrics:
            metric_func = METRICS.get(metric)
            if not metric_func:
                raise ValueError(f"Metric '{metric}' is not defined in METRICS.")
            # Check if additional parameters are needed based on the metric type
            if metric == "nselog_mm7q":
                e = metric_func(md, obs, dry_month, gauges, gid)
            elif metric in {"mae_peak_timing", "mape_peak_magnitude"}:
                e = metric_func(peaks[run], window)
            else:
                e = metric_func(md, obs, gauges, gid)
            # Special case for the 'kge' metric to extract specific component    
            if metric == "kge":
                e = e["kge"]
            # store the metric value
            run_results[metric] = e
        
        # store the results for each run
        all_run_results[run] = run_results
        print(f"Results for {run} calculated")
    
    # ORGANIZE THE RESULTS INTO A DATASET
    runs = list(all_run_results.keys())
    # Create a dictionary to store the data variables
    data_vars = {}
    # Populate the data variables
    for metric in metrics:
        data = np.array([all_run_results[run][metric] for run in runs])
        data_vars[metric] = (('runs', 'wflow_id'), data)
    # Create the xarray Dataset
    ds = xr.Dataset(
        data_vars=data_vars,
        coords={
            'runs': runs,
            'wflow_id': gauges
        },
        attrs={"metrics": metrics,
               "starttime": starttime,
               "endtime": endtime,
               "wflow_id": gauges}
    )
    
    ds.to_netcdf(Path(output_dir, f'{filename}.nc'))
    
    return ds


def normalize_mae(
    val:float,
    window: int,
    ):
    """
    window: Size of window to consider on each side of the observed peak for finding the simulated peak.
    Assuming less than 1 hour lead lag is as good as perfect
    """
    
    if val < 1:
        norm = 1
    
    #best values approach 1
    else:
        norm  = 1 - val / window
    
    return norm

def normalize_mape(val:float):
    """__summary__"""
    norm = 1 - val
    return norm

def weighted_euclidean(
    coef: tuple | list,
    weights: tuple | list,
    weighted = True,
):
    """_summary_"""

    if weighted and len(weights) != len(coef):
        raise ValueError("The length of weights should be equal to the length of coef")
    
    if weighted and sum(weights) != 1:
        raise ValueError("The sum of weights should be equal to 1")
    
    if not weighted:
        weights = [1] * len(coef)

    dist = [
        w * (1-item)**2 for item, w in zip(coef, weights) 
    ]

    res = np.sqrt(sum(dist))

    return np.array(res.round(4))


def calculate_weighted_euclidean(
    ds_fn: Path | str,
    weighted_metrics: tuple | list,
    weights: tuple | list,
    window: int,
):
    ds = xr.open_dataset(ds_fn)
    
    # Create a new DataArray to store euclidean values
    euclidean_values = xr.DataArray(
        np.zeros((len(ds.runs), len(ds.wflow_id))),
        coords=[('runs', ds.runs.values), ('wflow_id', ds.wflow_id.values)],
        dims=['runs', 'wflow_id']
    )
    
    # calculate the euclidean distance
    for run in ds.runs.values:
        for g in ds.wflow_id.values:
            coef = []
            for metric in weighted_metrics:
                value = ds.sel(runs=run, wflow_id=g)[metric].values
                if metric == "mae_peak_timing":
                    value = normalize_mae(value, window)
                elif metric == "mape_peak_magnitude":
                    value = normalize_mape(value)
                coef.append(value)
            
            euclidean_value = weighted_euclidean(np.array(coef), weights=weights, weighted=True)
            euclidean_values.loc[{'runs': run, 'wflow_id': g}] = euclidean_value
            
    # Add the new DataArray to the dataset
    ds['euclidean'] = euclidean_values
    
    return ds


if __name__ == "__main__":
    
    import platform
    if platform.system() == "Windows":
        DRIVE = "p:/"
        PLATFORM = "Windows"
    elif platform.system() == "Linux":
        DRIVE = "/p"
        PLATFORM = "Linux"
    
    work_dir = Path(DRIVE, '11209265-grade2023', 'wflow', 'RWSOS_Calibration', 'meuse_random_spider')
    
    GaugeToPlot = Path(work_dir,'data', '4-output','wflow_id_add_HBV_new.csv')

    starttime = '2005-08-01'
    endtime = '2018-02-22'
    output_dir = Path(work_dir, 'data/4-output').as_posix()
    
    ds_fn = work_dir / 'data/4-output/combined_output.nc'
    gauges = pd.read_csv(GaugeToPlot)['wflow_id'].tolist()
    metrics = ['kge', 'nse','nse_log','nselog_mm7q', 'mae_peak_timing', 'mape_peak_magnitude']
    dry_month = [6, 7, 8, 9, 10]
    window = 72
    
    # ds=final_performance(
    #     ds=ds_fn,
    #     gid='wflow_id',
    #     starttime=starttime,
    #     endtime=endtime,
    #     gauges=gauges,
    #     metrics=metrics,
    #     dry_month=dry_month,
    #     window=window,
    #     output_dir=output_dir,
    #     filename='final_performance_sobek_gauges',
    # )
    
    weighted_metrics = ["kge", "nselog_mm7q", "mae_peak_timing", "mape_peak_magnitude"]
    weights=[0.2, 0.25, 0.3, 0.25]
    
    ds_fn = Path(output_dir, 'final_performance_sobek_gauges.nc')
    
    ds = calculate_weighted_euclidean(
        ds_fn=ds_fn,
        weighted_metrics=weighted_metrics,
        weights=weights,
        window=72,
    )
    
    ds.to_netcdf(Path(output_dir, 'final_performance_sobek_gauges_add_euclidean.nc'))
    
    
    # CHOOSE THE TOP MODEL!
    # 1. BASED ON SOBEK GAUGES
    # Create a dataframe where rows are wflow_id, and columns are runs Top_1 to Top_10
    # and cell values are euclidean of that runs and wflow_id from ds
    df = pd.DataFrame(index=gauges, columns=['location']+ [f'Top_{i}' for i in range(1, 11)])
    df['location'] = pd.read_csv(GaugeToPlot)['location'].values
    for i in range(1, 11):
        run = f'Top_{i}'
        df[run] = ds.sel(runs=run).euclidean.values

    # check euclidean values
    # Calculate the average performance of each model
    average_performance = df.drop(columns='location').mean()
    # Count the number of locations where each model has the best (lowest) performance
    best_performance_count = (df.drop(columns='location').eq(df.drop(columns='location').min(axis=1), axis=0)).sum()

    print(average_performance, best_performance_count)
    
    # 2. BASED ON ALL GAUGES
    ds_fn_all = Path(output_dir, 'final_performance_allgauges.nc')
    ds_all = calculate_weighted_euclidean(
        ds_fn=ds_fn_all,
        weighted_metrics=weighted_metrics,
        weights=weights,
        window=72,
    )
    
    ds_all.to_netcdf(Path(output_dir, 'final_performance_allgauges_add_euclidean.nc'))
    
    df_all = pd.DataFrame(index=ds_all.wflow_id.values.tolist(), columns=[f'Top_{i}' for i in range(1, 11)])
    for i in range(1, 11):
        run = f'Top_{i}'
        df_all[run] = ds_all.sel(runs=run).euclidean.values
        
    average_performance_all = df_all.mean()
    best_performance_count_all = (df_all.eq(df_all.min(axis=1), axis=0)).sum()

    print(average_performance_all, best_performance_count_all)
    