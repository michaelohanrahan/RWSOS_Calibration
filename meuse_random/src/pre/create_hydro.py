from pathlib import Path

import numpy as np
import pyflwdir
import xarray as xr
from hydromt import DataCatalog, flw

def gen_hydro_from_dem(
    dem_cat: Path | str,
    output_dir: Path | str,
    catalog: Path | str,
    ds_name: str = "ned_hydro",
    upscale: int = None,
):
    """_summary_."""
    dc = DataCatalog(data_libs=[dem_cat])
    dem = dc.get_rasterdataset("ned")
    crs = dem.raster.crs
    epsg = crs.to_epsg()

    hydro_ds = dem.to_dataset(name="elevtn")
    hydro_ds["flwdir"] = flw.d8_from_dem(
        da_elv=dem,
        gdf_stream=None,
        max_depth=-1,  # no local pits
        outlets="edge",
        idxs_pit=None,
    )

    # Set the deminsions
    dims = hydro_ds.raster.dims

    # Create a PyFlwDir object from the dataset
    flwdir = flw.flwdir_from_da(hydro_ds["flwdir"])

    # uparea
    uparea = flwdir.upstream_area(unit="km2")
    hydro_ds["uparea"] = xr.Variable(dims, uparea, attrs=dict(_FillValue=-9999))

    if upscale is not None:
        ds_name += f"_ihu{upscale}"

        da, flwdir = flw.upscale_flwdir(
            hydro_ds,
            flwdir,
            upscale,
            method="ihu",
            uparea_name = "uparea",
        )

        hydro_ds = da.to_dataset().reset_coords(["x_out", "y_out", "idx_out"])
        hydro_ds.raster.set_crs(epsg)
        outidx = hydro_ds["idx_out"].values
        mask = np.where(dem.values==dem.raster.nodata)
        if len(mask[0]) > 0:
            novalidx = np.ravel_multi_index(
                ((mask[0][0],),(mask[1][0],)),
                dem.raster.shape
            )
            outidx = np.where(outidx==-1, novalidx, outidx)
        un_idx = np.unravel_index(outidx, dem.shape)
        hydro_ds["elevtn"] = xr.Variable(
            dims, 
            dem.values[un_idx[0], un_idx[1]],
            attrs=dict(_FillValue=dem.raster.nodata)
        )
        # Redo uparea
        uparea = flwdir.upstream_area(unit="km2")
        hydro_ds["uparea"] = xr.Variable(dims, uparea, attrs=dict(_FillValue=-9999))

        hydro_ds = hydro_ds.drop_vars(
            ["x_out", "y_out", "idx_out"]
        )

    # stream order
    strord = flwdir.stream_order()
    hydro_ds["strord"] = xr.Variable(dims, strord)
    hydro_ds["strord"].raster.set_nodata(255)

    # slope
    slope = pyflwdir.dem.slope(
        elevtn=hydro_ds["elevtn"].values,
        nodata=hydro_ds["elevtn"].raster.nodata,
        latlon=True,  # True if geographic crs, False if projected crs
        transform=hydro_ds["elevtn"].raster.transform,
    )
    hydro_ds["slope"] = xr.Variable(dims, slope)
    hydro_ds["slope"].raster.set_nodata(hydro_ds["elevtn"].raster.nodata)

    # basin at the pits locations
    basins = flwdir.basins(idxs=flwdir.idxs_pit).astype(np.int32)
    hydro_ds["basins"] = xr.Variable(dims, basins, attrs=dict(_FillValue=0))

    # basin index file
    gdf_basins = hydro_ds["basins"].raster.vectorize()

    hydro_ds

    # Exporting
    # Export the gridded data as tif files in a new folder
    # export the hydrography data as tif files (one per variable)
    hydro_ds.raster.to_mapstack(
        root=Path(output_dir, ds_name),
        driver="GTiff",
    )

    # export the basin index as geosjon
    gdf_basins.to_file(
        Path(output_dir, f"{ds_name}_basins.geojson"), driver="GeoJSON"
    )

    add_yml = f"""{ds_name}:
  crs: {epsg}
  data_type: RasterDataset
  driver: raster
  kwargs:
    chunks:
      x: 6000
      y: 6000
  meta:
    category: topography
    processing_notes: prepared from NED by deriving flow directions using pyflwdir.
    processing_script: create_hydro.py
  path: HYDRO/{ds_name}/*.tif
  rename:
    slope: lndslp

{ds_name}_index:
  crs: {epsg}
  data_type: GeoDataFrame
  driver: vector
  meta:
    category: topography
    processing_notes: prepared from NED by deriving flow directions using pyflwdir.
    processing_script: create_hydro.py
  path: HYDRO/{ds_name}_basins.geojson
  rename:
    value: basid"""

    with open(catalog, "a") as w:
        w.write("\n")
        w.write(add_yml)
        w.write("\n")

    pass

if __name__ == "__main__":
    gen_hydro_from_dem(
        "/p/1000365-002-wflow/tmp/usgs_wflow/data/DEM/dem.yml",
        "/p/1000365-002-wflow/tmp/usgs_wflow/data/HYDRO",
        "/p/1000365-002-wflow/tmp/usgs_wflow/data/catalog.yml",
    )