from pathlib import Path

import geopandas as gpd
import xarray as xr
from hydromt.raster import RasterDataset


def main(
    p: Path | str,
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    gauge_ids: tuple | list,
    sub_catch: Path | str,
    out: Path | str,
):
    """_summary_"""
    # Load the current main dataset
    ds = xr.open_dataset(p)

    # Make integers of the gauge_ids
    gauge_int = [
        int(item) for item in gauge_ids
    ]

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    # Get the relevant sub catchments
    vds = vds.astype({"value": int})
    vds = vds[vds.value.isin(gauge_int)]
    # Create the data mask
    mask = ds.raster.geometry_mask(vds)

    # Loop through the parameters
    for idx, (key, value) in enumerate(params.items()):
        vds[key] = value
        da = ds[params_lname[idx]]
        if params_method[idx] == "mult":
            da.values[mask] *= value
        elif params_method[idx] == "set":
            da.values[mask] = value          
        ds[params_lname[idx]] = da

    ds.to_netcdf(out)


if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]

        main(
            mod.params.dataset,
            mod.params.params,
            mod.params.params_lname,
            mod.params.params_method,
            mod.params.graph[mod.params.level]["elements"],
            mod.params.sub_catch,
            mod.output.staticmaps,
        )

    else:
        main(
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
            {"ksat": 5, "rd": 0.5, "st": 0.5},
            ["KsatHorFrac", "RootingDepth", "SoilThickness"],
            ["set", "mult", "mult"],
            ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000'],
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/out.nc"
        )
    pass


