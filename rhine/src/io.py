import zarr
from pathlib import Path

import xarray as xr


def export_to_zarr(
    ds: xr.Dataset,
    chunks: dict,
    out_fn: Path | str,
):
    """_summary_."""
    compressor = zarr.Blosc(cname="zstd", clevel=3, shuffle=1)

    encoding = {}
    ds.encoding = {}
    for v in ds.data_vars:
        ds[v].encoding = {}
        if compressor is not None:
            encoding[v] = {"compressor": compressor}

    ds.chunk(chunks).to_zarr(
        out_fn, 
        compute=True, 
        safe_chunks=True, 
        consolidated=False,
        encoding=encoding
    )

    zarr.convenience.consolidate_metadata(out_fn)


if __name__ == "__main__":
    ds = xr.open_mfdataset(
        "p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/HISTORY/cmip6_EcEarth/*.nc",
    )
    export_to_zarr(
        ds,
        "p:/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/HISTORY/cmip6_EcEarth.zarr",
    )
    pass