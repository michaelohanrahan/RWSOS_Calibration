from pathlib import Path

import geopandas as gpd
import xarray as xr
from hydromt.raster import RasterDataset
from setuplog import setup_logging
import shutil
import traceback


def main(
    l,
    p: Path | str,
    params: dict,
    params_lname: tuple | list,
    params_method: tuple | list,
    gauge_ids: tuple | list,
    sub_catch: Path | str,
    lakes_in: Path | str,
    lakes_out: Path | str,
    out: Path | str,
):
    """
    Load the dataset (staticmaps) then update params for the desired gauge_ids
    the paramslname should be the same as the dataset key
    #DONE: make sure the calibrecipe lnames match 
    #DONE: create Floodplain_N
    """
    # if not os.path.exists(out):
    ST_val = params["st"]
    l.info(f"parsed soil thickness value: {ST_val}")
    
    # Load the current main dataset
    ds = xr.open_dataset(p)
    l.info(f"Loading dataset from {p}")
    # Make integers of the gauge_ids
    gauge_int = [
        int(item) for item in gauge_ids
    ]
    l.info(f"Updating the following gauge_ids: {gauge_int}")

    # Load the geometries
    vds = gpd.read_file(sub_catch)
    
    # Get the relevant sub catchments
    #vds is now a geodataframe
    vds = vds.astype({"value": int})
    #vds is now a geodataframe with only the gauges of interest
    vds = vds[vds.value.isin(gauge_int)]
    
    # Create the data mask
    mask = ds.raster.geometry_mask(vds)

    # Loop through the parameters
    for idx, (key, value) in enumerate(params.items()):
        #The key is the var id
        vds[key] = value
        da = ds[params_lname[idx]]
        if params_method[idx] == "mult":
            da.values[mask] *= value
        elif params_method[idx] == "set":
            da.values[mask] = value
        elif params_method[idx] == "add":
            da.values[mask] += value          
        ds[params_lname[idx]] = da

    ds.to_netcdf(out)

    l.info(f"Writing dataset to {out}")

    #make sure the lakes h-q rating relation is also copied
    if not os.path.exists(lakes_out):
        if len(lakes_in)!=len(lakes_out):
            l.error(f"lakes_in and lakes_out should have the same length")
            raise ValueError("lakes_in and lakes_out should have the same length")
        for src, dst in zip(lakes_in, lakes_out):
            shutil.copy(src, dst)
            l.info(f"Copying {src} to {dst}")
        
        
    


if __name__ == "__main__":
    l=setup_logging('data/0-log', '03-set_calib_params.log')
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]

            main(
                l,
                mod.params.dataset,
                mod.params.params,
                mod.params.params_lname,
                mod.params.params_method,
                mod.params.graph[mod.params.level]["elements"],
                mod.params.sub_catch,
                mod.params.lakes_in,
                mod.output.lakes_hqs, 
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

    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e

