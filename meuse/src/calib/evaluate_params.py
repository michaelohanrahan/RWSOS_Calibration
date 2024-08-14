from pathlib import Path
import os
# import ast
import numpy as np
import pandas as pd
import xarray as xr
from setuplog import setup_logging
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
    random_df: Path | str,
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
    gauges = graph[level]["elements"]

    obs = xr.open_dataset(observed)
    l.info(f"Opened observed data from {observed}")
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
    
    l.info(f"Params: {params}")
    param_coords = create_coords_from_params(params)
    l.info(f"Created parameter coordinates")    
    l.info(f"Parameter coordinates: {param_coords}")
    
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
        upgauge_ids = graph_node[str(g)]['_deps']  # get the upstream gauges ids for gauge g
        
        for topx in _out_ds_sel.index:
            
            # search in random_df for the matching row for the Top_x
            _mask = pd.Series([True] * len(random_df))
            
            for key, value in _out_ds_sel[topx].items():
                _mask = _mask & (random_df[key] == value)
            random_df_sel = random_df[_mask]
            
            # select out the upstream gauges from random_df and fill in the _out_ds
            for upgauge in upgauge_ids:
                _out_ds.loc[(level, upgauge), topx] = random_df_sel[str(upgauge)].values[0]
                
    
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
    random_df = pd.read_csv(work_dir / 'random_df.csv', index_col=0)
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
        random_df=random_df,
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
    
    
    # if "snakemake" in globals():
    #     mod = globals()["snakemake"]
        
    #     main(
    #         mod.input,
    #         mod.params.observed_data,
    #         mod.params.graph[mod.params.level]["elements"],
    #         mod.params.params,
    #         mod.params.starttime,
    #         mod.params.endtime,
    #         mod.params.metrics,
    #         mod.params.weights,
    #         mod.output.best_params,
    #     )

    # else:
    #     from create_set_params import create_set
    #     from dependency_graph import sort_graph
    #     graph = sort_graph(
    #         Path(r"c:\git\puget\res\king_graph.json")
    #     )
    #     lnames, methods, ds = create_set(r"c:\git\puget\res\calib_recipe.json")

    #     main(
    #         modelled=[
    #             Path(
    #                 "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2", 
    #                 "ksat{}_rd{}_st{}".format(*item.values()), 
    #                 "run_default", 
    #                 "output_scalar.nc"
    #             ) 
    #             for item in 
    #             ds.to_dict(orient="records")
    #         ],
    #         observed="p:/1000365-002-wflow/tmp/usgs_wflow/data/GAUGES/discharge_obs_combined.nc",
    #         dry_month=None,
    #         window=None,
    #         gauges=graph["level2"]["elements"],
    #         params=ds.to_dict(orient="records"),
    #         starttime="2011-01-01T00:00:00",
    #         endtime="2011-12-31T00:00:00",
    #         metrics=["kge", "rld"],
    #         weights=[0.6, 0.4],
    #         out="p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level2/best_params.csv",
    #     )