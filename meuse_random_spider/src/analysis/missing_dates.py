from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


CLIM_SCENARIOS = [
    "cmip6_HadGemHH",
    "cmip6_HadGemHM",
    "cmip6_HadGemHMsst",
]


def main(
    p: Path | str,
):
    """_summary_."""
    p = Path(p)

    data = {}

    for clim in CLIM_SCENARIOS:
        ds = xr.open_zarr(
            Path(p, f"{clim}.zarr")
        )
        all_time = pd.date_range(
            start=ds.time.min().values, 
            end=ds.time.max().values, 
            freq="3h"
        )
        missing = all_time.difference(ds.time.values)

        missing = [
            item.strftime("%m-%dT%H:%M:%S") for item in missing
        ]
        missing = np.unique(missing)

        data[clim] = missing

    df = pd.DataFrame(data) 
    df.to_csv("./tmp/future_missing.csv")
    pass


if __name__ == "__main__":
    main(
        "p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/FUTURE",
    )
    ds = xr.open_zarr("p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/HISTORIC/cmip6_HadGemHH.zarr")
    pass