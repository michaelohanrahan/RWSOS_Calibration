from pathlib import Path

import xarray as xr
from hydromt import DataCatalog


def create_orography(
    da: xr.DataArray,
    ds_like: xr.Dataset,
    fname: str,
    out_dir: Path | str,
):
    """_summary_."""
    da = da.raster.reproject_like(
        ds_like,
        method="average",
    )

    da.name = fname
    ds = da.to_dataset()
    ds = ds.drop("height")
    ds = ds.isel(time=0).drop("time")
    ds.to_netcdf(
        Path(out_dir, "cmip6_GFDL_orog.nc"),
        encoding={fname: {"zlib": True}},
    )


def reorient_orography(
    p: Path | str,
    d: DataCatalog,
    fname: str,
):
    """_summary_."""
    p = Path(p)

    out_dir = p.parent

    for _f in p.iterdir():
        bname = _f.stem.rsplit("_",1)[0]
        ds = xr.open_dataset(_f)

        dname = f"{bname}_historic"
        if dname not in d:
            continue

        ds_like = d.get_rasterdataset(dname)    

        if ds.lon.values[0] > 180:
            ds = ds.assign_coords({"lon": ds.lon.values - 360})
        if ds.lat.values[-1] > ds.lat.values[0]:
            ds = ds.raster.flipud()

        ds = ds.raster.reproject_like(
            ds_like,
            method="nearest",
        )

        ds = ds.drop("height")
        ds = ds.rename_vars({"orog": fname})

        ds.to_netcdf(
            Path(out_dir, f"{bname}_orog.nc"),
            encoding={fname: {"zlib": True}}
        )
        pass
    pass


if __name__ == "__main__":
    # Read the data in
    d = DataCatalog(
        data_libs=[
            "deltares_data",
            "p:/1000365-002-wflow/tmp/usgs_wflow/data/catalog.yml",    
            "p:/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml"
        ]
    )
    cl = d.get_rasterdataset("cmip6_GFDL_historic")
    dem = d.get_rasterdataset(
        "era5_orography", 
        bbox=cl.raster.bounds, 
        buffer=2, 
    )

    # Reorient existing files
    # reorient_orography(
    #     "p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/OROG/RAW",
    #     d,
    #     "elevtn",
    # )

    # Create one for GFDL
    create_orography(
        dem,
        cl,
        "elevtn",
        "p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/OROG",
    )
    pass
