from pathlib import Path
import os
import numpy as np
import pandas as pd
import xarray as xr
from setuplog import setup_logging
import traceback
from timeit import default_timer as timer

from metrics import kge, nselog_mm7q, mae_peak_timing, mape_peak_magnitude, weighted_euclidean
from metrics import _obs_peaks, _sim_peaks


def create_index_from_params(params: dict) -> pd.Index:
    """
    Create a Pandas Index from a single set of parameters.
    
    Args:
        params (dict): Dictionary containing parameter names as keys and a single set of values.

    Returns:
        pd.Index: A Pandas Index object representing the single set of parameters.
    """
    keys = list(params.keys())
    params_values = list(params.values())
    
    # index = pd.MultiIndex.from_tuples([tuple(params_values)], names=keys)
    index = pd.Index([tuple(params_values)], name=tuple(keys))
    return index

def parse_params_from_path(file_path):
    # Extract the part of the path that contains the parameters
    path_parts = Path(file_path).parts
    params_dict = {part.split('~')[0]:np.float64(part.split('~')[1]) for part in path_parts if '~' in part}  # Assuming the parameters are in the third last part of the path
    return params_dict

def main(
    l,
    modelled: tuple | list,
    observed: Path | str,
    dry_month: list,
    window: int,
    gauges: tuple | list,
    params: dict | str,
    starttime: str,
    endtime: str,
    metrics: tuple | list,
    weights: tuple | list,
    out: Path | str,
    gid: str,
):
    """
    Perform evaluation of model parameters.

    Args:
        modelled (tuple | list): List of paths to modelled data files.
        observed (Path | str): Path to observed data file.
        dry_month (list): List of dry months.
        window (int): Window size for peak calculation.
        gauges (tuple | list): List of gauge IDs.
        params (dict | str): Parameters for the model.
        starttime (str): Start time for the evaluation period.
        endtime (str): End time for the evaluation period.
        metrics (tuple | list): List of metrics to evaluate.
        weights (tuple | list): List of weights for the metrics.
        out (Path | str): Output directory.
        gid (str): Gauge key for simulation data.
    """
    #original example indexed by 'Q_gauges_obs' and 'time'
    #we use gid 
    #our data is indexed by 'wflow_id' 'runs' and 'time'
    
    METRICS = {metric: globals()[metric] for metric in metrics}
    
    out_dir = Path(out).parent
    os.makedirs(out_dir, exist_ok=True)
    obs = xr.open_dataset(observed)
    obs = obs.sel(runs='Obs.', time=slice(starttime, endtime)) 

    metric_values = {
        item: [] for item in metrics
    }
    evals=[]

    _m=modelled 
    
    if not Path(_m).exists():
        l.warning(f"\n{'*'*10}\n{_m}\nDoes not exist \n{'*'*10}")
    
    # Open and slice temporally 
    try:
        md = xr.open_dataset(_m)
        
    except Exception as e:
        l.error(f"Failed to open modelled data from {_m}")
        l.error(traceback.format_exc())
    
    md = md.sel(time=slice(starttime, endtime))
    
    if len(md.time.values) == 0:
        l.warning(f"\n{'*'*10}\n{_m}\nIs not a complete time series, Skipping...\n{'*'*10}")
        with open(Path(out_dir, "failed.nc"), 'a') as f:
            f.write(f"")
        raise ValueError(f"{_m} is not a complete time series")
    
    params = parse_params_from_path(_m)
    
    if md[gid].dtype != int:
        md[gid] = md[gid].astype(np.int64)

    if any("peak" in metric for metric in metrics):
        start = timer()
        obs_peaks = {
            g: _obs_peaks(obs.sel(wflow_id=g).Q)
            for g in gauges
        }
        peaks = {
            g: _sim_peaks(sim=md.sel({gid: g}).Q, 
                          obs=obs_peaks[g],
                          window=window,)
            for g in gauges
        }
        end = timer()
        l.info(f"Calculated peaks in {end-start} seconds")
    else:
        peaks = None

    # Calculate the metric values
    for metric in metrics:
        start = timer()
        metric_func = METRICS.get(metric)
        if not metric_func:
            raise ValueError(f"Metric '{metric}' is not defined in METRICS.")
        
        # Check if additional parameters are needed based on the metric type
        if metric == "nselog_mm7q":
            e = metric_func(md, obs, dry_month, gauges, gid)
        elif metric in {"mae_peak_timing", "mape_peak_magnitude"}:
            e = metric_func(peaks, window)
        else:
            e = metric_func(md, obs, gauges, gid)
        
        # Special case for the 'kge' metric to extract specific component    
        if metric == "kge":
            e = e["kge"]
        
        metric_values[metric].append(e)
        end = timer()
        l.info(f"Calculated {metric} in {end-start} seconds")
        l.info(f"{metric} value: {e}")
        evals.append(e)
    
    l.info(f"Calculated metrics for {_m}")
    
    # Get the euclidean distance
    res = weighted_euclidean(
        np.array(evals),
        weights=weights,
        weighted=True,
    )

    param_coords = create_index_from_params(params)
    
    ds = None
    for metric in metrics:
        l.info(metric)
        #metric_values[metric] has an array of len 60 at level 0
        da = xr.DataArray(
            metric_values[metric], 
            coords=[('params', param_coords), ('gauges', gauges)], 
            dims=['params', 'gauges'],
            attrs={"metric": metric}
        )
        l.info(da)
        da.name = metric
        if ds is None:
            ds = da.to_dataset()
            continue

        ds[metric] = da
    
    # l.info(ds)
    
    if res.ndim == 1:
        res = res.reshape((len(param_coords), len(gauges)))
    
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

    ds = ds.unstack()
    ds.to_netcdf(
        Path(out)
    )
    l.info(f"Saved the performance dataset to {out_dir}")
    #
    
    with open(out_dir / "evaluate.done", "w") as f:
        f.write("")
    


if __name__ == "__main__":
    
    """
    This module is used to evaluate parameters for a model. 

    it will evaluate per run returning a netcdf file with the performance metrics for each run.

    """
    
    l = setup_logging("data/0-log", "04-evaluate_param_per_run.log")
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            main(
                l,
                modelled=mod.input.sim,
                observed=mod.params.observed_data,
                dry_month=mod.params.dry_month,
                window=mod.params.window,
                gauges=mod.params.graph[f"level{mod.params.level}"]["elements"],
                params=mod.params.params,
                starttime=mod.params.starttime,
                endtime=mod.params.endtime,
                metrics=mod.params.metrics,
                weights=mod.params.weights,
                out=mod.output.performance,
                gid=mod.params.gaugeset,
            )

        else:
            from create_set_params import create_set
            import json
            import yaml            
            from numpy import random
            
            gpath = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\Hall_levels_graph.json"
            
            with open(gpath) as f:
                graph = json.load(f)
            
            cfg = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\config\calib.yml"
            
            with open(cfg) as f:
                cfg = yaml.safe_load(f)
            
            lnames, methods, ds = create_set(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\config\MINIMAL_calib_recipe.json")
            
            #random integer
            randint = random.randint(0, len(ds))
            
            #random param combo
            params = "".join([f"{col}~{val}/" for col, val in zip(ds.columns, ds.iloc[randint].values)])
            
            modelled = rf"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/{params}/output_scalar.nc"
            print("Modelled file exists: ", Path(modelled).exists())
            print("Metrics: ", cfg["metrics"][0:2])
            print("Weights: ", cfg["weights"][0:2])
            main(l,
                modelled=modelled,
                observed=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external\discharge_hourlyobs_HBV_combined.nc",
                dry_month=[6, 7, 8],
                window=60,
                gid='Q_gauges_Hall',
                gauges=graph["level0"]["elements"],
                params=params,
                starttime= "2008-01-01T01:00:00.000000000",
                endtime="2018-02-22T00:00:00.000000000", 
                metrics=cfg["metrics"][0:2],
                weights=[0.5,0.5],
                out=rf"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/{params}/evaluated.nc",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise e
    