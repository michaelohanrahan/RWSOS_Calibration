from pathlib import Path

import numpy as np
import xarray as xr
from hydromt import DataCatalog

from puget import CMIP6_MODELS, export_to_zarr

ERA5_DROPVARS = ['cape', 'temp_dew', 'tcwv', 'kout', 'wind10_u', 'wind10_v']
ERA5_VAR_METHOD = {
    "precip": np.sum,
    "press_msl": np.mean,
    "kin": np.mean,
    "temp": np.mean,
}

def main(
    dc: DataCatalog,
    cl_name: str,
    stime: str,
    etime: str,
    output_dir: Path | str, 
):
    """_summary_."""

    # Setup the data
    ds_like = dc.get_rasterdataset(cl_name)
    ds = dc.get_rasterdataset(
        "era5_hourly_zarr", 
        bbox=ds_like.raster.bounds, 
        buffer=2, 
        time_tuple=(stime, etime)
    )
    ds = ds.drop_vars(ERA5_DROPVARS)
    ds = ds.raster.reproject_like(
        ds_like
    )

    # Set an empty outgoing dataset
    out_ds = None
    out_name = cl_name.rsplit("_",1)[0]
    for var in ds.data_vars:
        da = ERA5_VAR_METHOD[var](
            ds[var].resample(
                time="3H",
                closed="right",
                label="right",
            )
        )
        if out_ds is None:
            out_ds =  da.to_dataset()
            continue
        out_ds[var] = da
        pass
    out_ds = out_ds.assign_attrs(ds.attrs)
    export_to_zarr(
        out_ds,
        chunks={"time":3000},
        out_fn=Path(output_dir, f"{out_name}.zarr")
    )


if __name__ == "__main__":
    # get the datacatalog and datasets
    dc = DataCatalog(
        data_libs=[
            "/p/1000365-002-wflow/tmp/usgs_wflow/data/deltares_data.yml",
            "/p/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml",
        ]
    )
    starttime = "1969-12-31T22:00:00"
    endtime = "2014-12-31T21:00:00"
    # execute the main function
    for cl_name in CMIP6_MODELS:     
        main(
            dc,
            cl_name=f"{cl_name}_orog",
            stime=starttime,
            etime=endtime,
            output_dir="/p/1000365-002-wflow/tmp/usgs_wflow/data/TMP/era5_climate"
        )