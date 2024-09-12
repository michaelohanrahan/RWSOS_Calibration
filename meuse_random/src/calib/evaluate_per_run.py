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
from filelock import FileLock


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
    modelled: Path | str,
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
        l.error(e)
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
    
    md, obs = xr.align(md,obs, join='inner')
    
    l.info(f"sim time: {md.time}")
    l.info(f"obs time: {obs.time}")
    
    
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
    
    # Extract the level from the modelled file path
    # allowing easy access to eucldean 
    level = None
    for part in Path(_m).parts:
        if "level" in part:
            level = part
            break
    split = str(out).split(level)[0]
    results_file = Path(split) / level / f"results_{level}.txt"
    l.info(f"appending results to: {results_file}")
     
    append_results_to_file(results_file, gauges, _m, res, l)  # Call the new function

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


def append_results_to_file(results_file: Path, gauges: list, _m: Path, res: np.ndarray, l) -> None:
    """Append results to the specified results file with file locking."""
    # Use FileLock to ensure safe file access
    lock_path = results_file.with_suffix('.lock')
    with FileLock(str(lock_path)):
        # Append results to the specified results file
        with open(results_file, 'a') as f:
            # Write the header if the file is empty
            if os.stat(results_file).st_size == 0:
                header = "params," + ",".join(map(str, gauges)) + "\n"
                f.write(header)

            # Write the file path and distances
            f.write(f"{parse_params_from_path(Path(_m))}," + ",".join(map(str, res)) + "\n")


if __name__ == "__main__":
    
    """
    This module is used to evaluate parameters for a model. 

    it will evaluate per run returning a netcdf file with the performance metrics for each run.
    
    This can be grouped across multiple cpus and thus is much more efficient than the old file loop. 

    """
    
    l = setup_logging("data/0-log", "04-evaluate_param_per_run.log")
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            main(
                l,
                modelled=mod.input.sim,
                observed=mod.params.observed,
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
            import json
            import yaml            
            gpath = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/Hall_levels_graph.json"
            
            with open(gpath) as f:
                graph = json.load(f)
            
            cfg = "/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse_random/config/calib.yml"
            
            with open(cfg) as f:
                cfg = yaml.safe_load(f)
            
            params = "ksat~1.06/f~1.5/rd~0.39/st~0.77/nr~1.0/ml~0.47/nl~1.02/nf~1.02"

            params_d = {str(p.split("~")[0]):float(p.split("~")[1]) for p in params.split('/')}

            modelled = f"/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/{params}/output_scalar.nc"

            main(l,
                modelled=modelled,
                observed="/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/1-external/discharge_hourlyobs_smoothed.nc",
                dry_month=[6, 7, 8],
                window=60,
                gid='Q_gauges_Hall',
                gauges=graph["level0"]["elements"],
                params=params_d,
                starttime= cfg["starttime"],
                endtime=cfg["endtime"], 
                metrics=cfg["metrics"],
                weights=cfg["weights"],
                out="/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/ksat~1.06/f~1.5/rd~0.39/st~0.77/nr~1.0/ml~0.47/nl~1.02/nf~1.02/evaluated.nc",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise e
