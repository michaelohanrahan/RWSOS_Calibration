from pathlib import Path

import pyflwdir
import xarray as xr
from hydromt import open_raster

def calc(
    p: Path | str,
):
    """_summary_."""

    r = open_raster(p)
    r.name = "elevtn"
    hydro_ds = r.to_dataset()
    slope = pyflwdir.dem.slope(
        elevtn=r.values,
        nodata=r.raster.nodata,
        latlon=True,  # True if geographic crs, False if projected crs
        transform=r.raster.transform,
    )
    hydro_ds["slope"] = xr.Variable(hydro_ds.raster.dims, slope)
    hydro_ds["slope"].raster.set_nodata(hydro_ds["elevtn"].raster.nodata)
    slope = hydro_ds["slope"]
    pass

if __name__ == "__main__":
    p = r"p:\1000365-002-wflow\tmp\usgs_wflow\data\HYDRO\ned_hydro_ihu8\elevtn.tif"
    calc(p)