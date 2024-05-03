from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from metrics import kge, rld, peakdis, weighted_euclidean


METRICS = {
    "kge": kge,
    "rld": rld,
    "pds": peakdis,
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
    gauges: tuple | list,
    params: tuple | list,
    starttime: str,
    endtime: str,
    metrics: tuple | list,
    weights: tuple | list,
    out: Path | str,
):
    """_summary_"""
    # output directory
    out_dir = Path(out).parent

    # Open the observed data and translate to m3/s
    obs = xr.open_dataset(observed)
    obs.Q.values = obs.Q.values * (0.3048**3)
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
            e = METRICS[metric](
                md,
                obs,
                gauges
            )
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

    # TODO create overal performance netcdf file
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

    # Best parameters sets by indexing minimum euc. dist.
    best_params = np.argmin(res, axis=0)

    _out = []
    for idx, best in enumerate(best_params):
        _out.append(params[best])
    
    # Create a pandas dataframe for the best parameters
    out_ds = pd.DataFrame(
        _out,
        index=gauges,
    )

    # Set the index name of the csv
    out_ds.index.name = "gauges"
    # Write to file
    out_ds.to_csv(out)


if __name__ == "__main__":
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