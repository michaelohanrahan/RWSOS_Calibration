from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

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

    obs = xr.open_dataset(observed)
    obs = obs.sel(time=slice(starttime, endtime)) 

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

        # Calculate the metric values
        for metric in metrics:
            metric_func = METRICS.get(metric)
            if not metric_func:
                raise ValueError(f"Metric '{metric}' is not defined in METRICS.")
            
            # Check if additional parameters are needed based on the metric type
            if metric == "nselog_mm7q":
                e = metric_func(md, obs, dry_month, gauges)
            elif metric in {"mae_peak_timing", "mape_peak_magnitude"}:
                e = metric_func(md, obs, window, gauges)
            else:
                e = metric_func(md, obs, gauges)
            
            # Special case for the 'kge' metric to extract specific component    
            if metric == "kge":
                e = e["kge"]
            
            metric_values[metric].append(e)
            evals.append(e)

        # Get the euclidean distance
        r = weighted_euclidean(
            np.array(evals),
            weights=weights,
            weighted=True,
        )
        res.append(r)

    param_coords = create_coords_from_params(params)
    ds = None
    for metric in metrics:
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
    ds.to_netcdf(
        Path(out_dir, "performance.nc")
    )

    #TODO: 
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
    print(out)
    # out_ds.to_csv(out)


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
    
    
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        main(
            mod.input,
            mod.params.observed_data,
            mod.params.graph[mod.params.level]["elements"],
            mod.params.params,
            mod.params.starttime,
            mod.params.endtime,
            mod.params.metrics,
            mod.params.weights,
            mod.output.best_params,
        )

    else:
        from create_set_params import create_set
        from dependency_graph import sort_graph
        graph = sort_graph(
            Path(r"c:\git\puget\res\king_graph.json")
        )
        lnames, methods, ds = create_set(r"c:\git\puget\res\calib_recipe.json")

        main(
            modelled=[
                Path(
                    "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2", 
                    "ksat{}_rd{}_st{}".format(*item.values()), 
                    "run_default", 
                    "output_scalar.nc"
                ) 
                for item in 
                ds.to_dict(orient="records")
            ],
            observed="p:/1000365-002-wflow/tmp/usgs_wflow/data/GAUGES/discharge_obs_combined.nc",
            dry_month=None,
            window=None,
            gauges=graph["level2"]["elements"],
            params=ds.to_dict(orient="records"),
            starttime="2011-01-01T00:00:00",
            endtime="2011-12-31T00:00:00",
            metrics=["kge", "rld"],
            weights=[0.6, 0.4],
            out="p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2/best_params.csv",
        )