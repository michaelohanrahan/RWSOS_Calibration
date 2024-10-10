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

def parse_params_from_path(file_path):
    # Extract the part of the path that contains the parameters
    path_parts = Path(file_path).parts
    params_dict = {part.split('~')[0]:np.float64(part.split('~')[1]) for part in path_parts if '~' in part}  # Assuming the parameters are in the third last part of the path
    return params_dict

def main(
    l,
    l,
    modelled: tuple | list,
    observed: Path | str,
    random_params: Path | str,
    best_10params_previous: Path | str,	
    dry_month: list,
    window: int,
    level: str,
    graph: dict,
    graph_node: dict,
    params: tuple | list,
    starttime: str,
    endtime: str,
    metrics: tuple | list,
    weights: tuple | list,
    gid: str,
    out: Path | str,
    gid: str,
):
    """
    Perform evaluation of model parameters.

    Args:
        modelled (tuple | list): List of paths to modelled data files.
        observed (Path | str): Path to the observed data file.
        level (str): Level of the graph.
        graph (dict): Graph of the model.
        # gauges (tuple | list): List of gauge names.
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
    
    # laod random_params
    random_params = pd.read_csv(random_params, index_col=0)

    res = []

    metric_values = {
        item: [] for item in metrics
    }

    # Create a list to store the parameters, using only non-empty runs to do so
    params = []
    path_ex = 0
    val_err = 0
    no_time = 0
    #failed runs
    
    # Loop through the modelled data
    for _m in modelled:
        evals = []
        
        if not Path(_m).exists():
            l.warning(f"\n{'*'*10}\n{_m}\nDoes not exist, Skipping...\n{'*'*10}")
            path_ex += 1
            continue
        
        # Open and slice temporally 
        try:
            md = xr.open_dataset(_m)
        except Exception as e:
            l.error(f"Failed to open modelled data from {_m}")
            l.error(traceback.format_exc())
            val_err += 1
            continue
        
        md = md.sel(time=slice(starttime, endtime))
        
        if len(md.time.values) == 0:
            l.warning(f"\n{'*'*10}\n{_m}\nIs not a complete time series, Skipping...\n{'*'*10}")
            no_time += 1
            continue
        
        params.append(parse_params_from_path(_m))
        
        if md[gid].dtype != int:
            md[gid] = md[gid].astype(np.int64)

        # Calculate the metric values
        for metric in metrics:
            metric_func = METRICS.get(metric)
            if not metric_func:
                raise ValueError(f"Metric '{metric}' is not defined in METRICS.")
            
            # Check if additional parameters are needed based on the metric type
            if metric == "nselog_mm7q":
                e = metric_func(md, obs, dry_month, gauges, gid, gid)
            elif metric in {"mae_peak_timing", "mape_peak_magnitude"}:
                e = metric_func(md, obs, window, gauges, gid, gid)
            else:
                e = metric_func(md, obs, gauges, gid, gid)
            
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

    #loop summary
    l.info(f"Failed to open {path_ex + val_err + no_time} modelled data files")
    l.info(f"- {path_ex} files did not exist")
    l.info(f"- {val_err} files failed to open, (probably 'nc4 error', file not valid)")
    l.info(f"- {no_time} files did not have a valid time dimension (empty)")
    l.info(f"Successfully opened and evaluated {len(params)} modelled data files")
    
    param_coords = create_coords_from_params(params)
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
    l.info(f"Created xarray dataset with metrics and euclidean distance")
    l.info(f"{ds}")
    ds.to_netcdf(
        Path(out_dir, "performance.nc")
    )
    l.info(f"Saved the performance dataset to {out_dir}")

    best_10params = np.argsort(res, axis=0)[:10]

    _out = []
    for idx, best in enumerate(best_10params.T):
        _out.append([params[i] for i in best])
    
    _out_ds = pd.DataFrame(
            _out,
            index=pd.MultiIndex.from_product([[level], gauges], names=['level', 'gauges']),
            columns=[f"Top_{i+1}" for i in range(10)]
        )
    
    # add for upstream gauges
    for g in gauges:
        _out_ds_sel = _out_ds.loc[(level, g)]  # select the top 10 parameters for the gauge g
        upgauge_ids = graph_node[str(g)]['_deps']  # get all the upstream gauges (not only direct upstream) ids for gauge g
        
        for topx in _out_ds_sel.index:
            
            # search in random_params for the matching row for the Top_x
            _mask = pd.Series([True] * len(random_params))
            
            for key, value in _out_ds_sel[topx].items():
                _mask = _mask & (random_params[key] == value)
            random_params_sel = random_params[_mask]
            
            # select out the upstream gauges from random_params and fill in the _out_ds
            for upgauge in upgauge_ids:
                _out_ds.loc[(level, upgauge), topx] = random_params_sel[str(upgauge)].values[0]
                
    
    # save the best 10 parameters result to dataframe (Load the previous level best parameters dataframe if exists)
    if level == 'level0' or not os.path.exists(best_10params_previous):
        out_ds = _out_ds
    else:
        prev_best_10params = pd.read_csv(best_10params_previous, index_col=[0, 1])
        out_ds_temp = prev_best_10params.copy()
        out_ds = pd.concat([out_ds_temp, _out_ds], axis=0)
    
    out_ds.to_csv(out)
    
    return ds, out_ds


if __name__ == "__main__":
    
    #NOTE: the TEST DATA cannot properly execute the following test code

    work_dir = Path(r"c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA")
    random_params = pd.read_csv(work_dir / 'random_df.csv', index_col=0)
    best_10params_previous = work_dir / 'best_10params_level4.csv'
    
    # import necessary variables
    import pickle as pk
    with open(work_dir/'create_set_params.pkl', 'rb') as f:
        dict = pk.load(f)
    lnames = dict['lnames']
    methods = dict['methods']
    df = dict['ds']
    
    # create inputs  
    modelled = [work_dir/'output_scaler_test1.nc', work_dir/'output_scaler_test2.nc']*5 # to generate best 1o param sets
    observed = work_dir / 'observed_data.nc'
    dry_month = [6, 7, 8, 9, 10]
    window = 72
    level = 'level5'
    import json
    graph = json.load(open(Path(work_dir / 'Hall_levels_graph.json')))
    graph_node = json.load(open(Path(work_dir / 'Hall_nodes_graph.json')))
    params = df.to_dict(orient="records")[:10]  # only use the first 10 records, because we only have 10 modelled results to test
    starttime = '2015-01-01T02:00:00'
    endtime = '2018-02-21T23:00:00'
    metrics = ["kge", "nselog_mm7q", "mae_peak_timing", "mape_peak_magnitude"]
    weights = [0.2, 0.25, 0.3, 0.25]
    out = work_dir / f'best_10params.csv' 
    
    # call function main()
    l = setup_logging(work_dir, "evaluate_params_0813.log")
    ds_performance, df_best_10params = main(
        l,
        modelled=modelled,
        observed=observed,
        random_params=random_params,
        best_10params_previous=best_10params_previous,
        dry_month=dry_month,
        window=window,
        level=level,
        graph=graph,
        graph_node=graph_node,
        params=params,
        starttime=starttime,
        endtime=endtime,
        metrics=metrics,
        weights=weights,
        gid='wflow_id',
        out=out,
    )
    
    
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
            Path("c:/CODING/NONPRODUCT/puget/res/king_graph.json")
        )
        lnames, methods, ds = create_set("c:/CODING/NONPRODUCT/puget/res/calib_recipe.json")

        main(
            [
                Path(
                    "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2", 
                    "ksat{}_rd{}_st{}".format(*item.values()), 
                    "run_default", 
                    "output_scalar.nc"
                ) 
                for item in 
                ds.to_dict(orient="records")
            ],
            "p:/1000365-002-wflow/tmp/usgs_wflow/data/GAUGES/discharge_obs_combined.nc",
            graph["level2"]["elements"],
            ds.to_dict(orient="records"),
            "2011-01-01T00:00:00",
            "2011-12-31T00:00:00",
            ["kge", "rld"],
            [0.6, 0.4],
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2/best_params.csv",
        )