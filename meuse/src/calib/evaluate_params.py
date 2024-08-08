from pathlib import Path
import os
import numpy as np
import pandas as pd
import xarray as xr
from setuplog import setup_logging
import traceback

#TODO: add peak timing to metrics, conforming to data structure
from metrics import kge, nselog_mm7q, mae_peak_timing, mape_peak_magnitude, weighted_euclidean

# define a dictionary to specify the metrics to use
METRICS = {
    "kge": kge,
    "nselog_mm7q": nselog_mm7q,
    "mae_peak_timing": mae_peak_timing,
    "mape_peak_magnitude": mape_peak_magnitude,
} 


def create_coords_from_params(
    params: dict,
):
    """_summary_"""
    keys = list(params[0].keys())

    params_values = list(
        zip(*[_dict.values() for _dict in params]),
    )

    coord = pd.MultiIndex.from_arrays(
        params_values, 
        names=keys
    )

    return coord


def main(
    l,
    modelled: tuple | list,
    observed: Path | str,
    dry_month: list,
    window: int,
    gauges: tuple | list,
    params: tuple | list,
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
        observed (Path | str): Path to the observed data file.
        gauges (tuple | list): List of gauge names.
        params (tuple | list): List of parameter values.
        starttime (str): Start time for data selection.
        endtime (str): End time for data selection.
        metrics (tuple | list): List of metrics to calculate.
        weights (tuple | list): List of weights for weighted metrics.
        out (Path | str): Path to the output file.

    Returns:
        None
    """
    #original example indexed by 'Q_gauges_obs' and 'time'
    #our data is indexed by 'wflow_id' 'runs' and 'time'
    #TODO: optionality to have additional index for alternative dimension 
    # output directory
    out_dir = Path(out).parent
    os.makedirs(out_dir, exist_ok=True)
    obs = xr.open_dataset(observed)
    l.info(f"Opened observed data from {observed}")
    obs = obs.sel(runs='Obs.', time=slice(starttime, endtime)) 

    res = []

    metric_values = {
        item: [] for item in metrics
    }

    # Loop through the modelled data
    for _m in modelled:
        evals = []

        # Open and slice temporally 
        md = xr.open_dataset(_m)
        md = md.sel(time=slice(starttime, endtime))
        
        if len(md.time.values) == 0:
            l.warning(f"\n{'*'*10}\n{_m}\nIs not a complete time series, Skipping...\n{'*'*10}")
            continue
        
        if md[gid].dtype != int:
            md[gid] = md[gid].astype(np.int64)

        l.info(f"Opened modelled data from {_m}")
        # Calculate the metric values
        for metric in metrics:
            metric_func = METRICS.get(metric)
            if not metric_func:
                raise ValueError(f"Metric '{metric}' is not defined in METRICS.")
            
            # Check if additional parameters are needed based on the metric type
            if metric == "nselog_mm7q":
                e = metric_func(md, obs, dry_month, gauges, gid)
            elif metric in {"mae_peak_timing", "mape_peak_magnitude"}:
                e = metric_func(md, obs, window, gauges, gid)
            else:
                e = metric_func(md, obs, gauges, gid)
            
            # Special case for the 'kge' metric to extract specific component    
            if metric == "kge":
                e = e["kge"]
            
            metric_values[metric].append(e)
            evals.append(e)
        l.info(f"Calculated metrics for {_m}")
        l.info(f"Metrics: {evals}")

        # Get the euclidean distance
        r = weighted_euclidean(
            np.array(evals),
            weights=weights,
            weighted=True,
        )
        l.info(f"Euclidean distance: {r}")
        res.append(r)
        l.info(f"\n{'='*30}\n -- Delete line 134 in eval when testing is complete -- \n{'='*30}\n")
        #getout
        break
    
    l.info(f"Params: {params}")
    param_coords = create_coords_from_params(params)
    l.info(f"Created parameter coordinates")    
    l.info(f"Parameter coordinates: {param_coords}")
    
    ds = None
    for metric in metrics:
        #mmetric_values[metric] has an array of len 60
        #there are now 128 values in the potential params
        da = xr.DataArray(
            metric_values[metric], 
            coords=[('params', param_coords), ('gauges', gauges)], 
            dims=['params', 'gauges']
        )
        da.name = metric
        if ds is None:
            ds = da.to_dataset()
            continue

        ds[metric] = da

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
    l.info(f"Created xarray dataset with metrics and euclidean distance")
    l.info(f"{ds}")
    ds.to_netcdf(
        Path(out_dir, "performance.nc")
    )
    l.info(f"Saved the performance dataset to {out_dir}")

    #TODO: test this best 10 params part
    # Best 10 parameters sets by indexing minimum euc. dist.
    best_10params = np.argsort(res, axis=0)[:10]

    _out = []
    for idx, best in enumerate(best_10params.T):
        _out.append([params[i] for i in best])
    
    # Create a pandas dataframe for the best parameters
    out_ds = pd.DataFrame(
        _out,
        index=gauges,
        columns=[f"Top_{i+1}" for i in range(10)]
    )

    # Set the index name of the csv
    out_ds.index.name = "gauges"
    # Write to file
    out_ds.to_csv(out)


if __name__ == "__main__":
    
    """
    This module is used to evaluate parameters for a model. 

    It contains a main function that takes in a variety of inputs including observed data, 
    graph elements, parameters, start and end times, metrics, weights, and a path to output 
    the best parameters. 

    The main function is designed to be run either directly or through Snakemake. 
    If run directly, it will use default paths and parameters defined in the script. 
    If run through Snakemake, it will use the input, params, and output defined in the 
    Snakemake rule.

    Functions:
    ----------
    main : Function to evaluate parameters and output the best ones.

    Modules:
    --------
    create_set_params : Module to create a set of parameters.
    dependency_graph : Module to sort the graph.

    """
    
    l = setup_logging("data/0-log", "04-evaluate_params.log")
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            main(
                l,
                modelled=mod.input,
                observed=mod.params.observed_data,
                dry_month=mod.params.dry_month,
                window=mod.params.window,
                gauges=mod.params.graph[f"level{mod.params.level}"]["elements"],
                params=mod.params.params,
                starttime=mod.params.starttime,
                endtime=mod.params.endtime,
                metrics=mod.params.metrics,
                weights=mod.params.weights,
                out=mod.output.best_10params,
                gid=mod.params.gaugeset,
            )

        else:
            from create_set_params import create_set
            from dependency_graph import sort_graph
            import json
            
            gpath = r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\Hall_levels_graph.json"
            with open(r) as f:
                graph = json.load(f)
            
            lnames, methods, ds = create_set(r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\config\TEST_calib_recipe.json")

            main(l,
                modelled=['p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~0.8/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~0.8/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~1.5/nr~0.5/rd~1.2/ml~0.6/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~0.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~0.5/st~1.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~2.0/f~1.5/st~0.5/nr~0.5/rd~0.8/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~0.5/nr~0.5/rd~1.2/ml~0.0/nl~1.2/nf~1.2/output_scalar.nc', 'p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.6/nl~0.8/nf~0.8/output_scalar.nc'],
                observed=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external\discharge_hourlyobs_HBV_combined.nc",
                gauges=graph["level0"]["elements"],
                params=ds.to_dict(orient="records"),
                starttime="2005-01-01T00:00:00",
                endtime="2006-12-31T00:00:00",
                metrics=["kge", "rld"],
                weights=[0.6, 0.4],
                out="p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2/best_params.csv",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise e
    