import xarray as xr
from pathlib import Path
from setuplog import setup_logging
import traceback
import pandas as pd


def main(
    l,
    files: list,
    out_nc: str|Path,
    best_n: int = 10,
):
    l.info(f"Combining {len(files)} files to {out_nc}")
    ds = xr.open_mfdataset(files, combine="by_coords")
    ds.to_netcdf(Path(out_nc))
    
    # Automatically detect parameter dimensions by excluding known non-parameter dimensions
    non_param_dims = ['gauges']
    param_dims = [dim for dim in ds.dims if dim not in non_param_dims]
    # print(param_dims)
    
    # Flatten the parameter dimensions into a single dimension
    flattened = ds.stack(params=param_dims).compute()
    print(flattened)
    # Rank the 'euclidean' values per gauge
    ranked = flattened['euclidean'].rank(dim='params')

    best_params = []
    for gauge in ds['gauges'].values:
        # Drop NaN values for this specific gauge to retain only valid parameter sets
        valid_euclidean = flattened['euclidean'].sel(gauges=gauge).dropna(dim='params')
        valid_ranked = ranked.sel(gauges=gauge).dropna(dim='params')
        top_n = min(best_n, len(valid_euclidean))
        
        best_for_gauge = {'gauge': gauge}
        for idx in range(top_n):
            param_vals = valid_ranked.where(valid_ranked == idx + 1, drop=True).params
            param_dict = {param: param_vals[param].values[0] for param in param_dims}
            best_for_gauge[f'Top_{idx+1}'] = param_dict
        best_params.append(best_for_gauge)
    
    
    final_df = pd.DataFrame(best_params)
    
    folder = Path(out_nc).parent
    out = Path(folder, "best_params").with_suffix('.csv')
    final_df.to_csv(out, index=False)
    # l.info(f"Processed and combined data saved to {out}")
    
    with open(Path(out_nc).parent / "eval.done", "w") as f:
        f.write("")
    l.info(f"Saved {out_nc}")


if __name__ == "__main__":
    l = setup_logging("data/0-log", "05-combine_evaluated.log")
    try:
        if 'snakemake' in globals():
            main(
                l,
                snakemake.input.performance_files,
                snakemake.output.performance,
            )
        else:
            main(l,
                [
                    r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\ksat~0.5\f~0.75\rd~1.0\st~0.5\nr~0.5\ml~0.3\evaluated.nc",
                    r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\ksat~0.75\f~0.5\rd~1.75\st~1.0\nr~1.0\ml~0.3\evaluated.nc",
                ],
                out_nc=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\performance_test.nc",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise