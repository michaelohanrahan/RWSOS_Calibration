from pathlib import Path

import geopandas as gpd
import xarray as xr
from hydromt.raster import RasterDataset
from setuplog import setup_logging
import shutil
import time
import traceback
import os


def main(
    l,
    base_grid: Path | str,
    level: int,
    ST_str: float,
    ST_key: str, 
    gauge_ids: list | tuple,
    sub_catch: Path | str,
    lakes_in: list | tuple,
):
    """
    """
    base_dir = base_grid.parent
    base_dir = base_dir.parent
    
    if len(ST_str) == 2:
        _ST = float(f"{ST_str[0]}.{ST_str[1]}")
    else:
        raise ValueError("ST_str must be a string of length 2")
        
    st_str = str(_ST).replace('.', '')
    
    # Create the output path
    sml = Path(base_dir, 'instates', f"staticmaps_L{level}_ST{st_str}.nc")
    
    # Copy the base_grid to a temporary file
    temp_base_grid = Path(base_dir, 'instates', f"temp_staticmaps_L{level}_ST{st_str}.nc")
    if not sml.exists():
        if not temp_base_grid.exists():
            shutil.copy(base_grid, temp_base_grid)
        
        if not temp_base_grid.exists():
            raise FileNotFoundError(f"Failed to copy {base_grid} to {temp_base_grid}")
        
        l.info(f"Copied {base_grid} to {temp_base_grid}")
        
        with xr.open_dataset(temp_base_grid) as _r:
            ds = _r.load()
        
        l.info(f"Loaded dataset from {temp_base_grid}")
        
        # Make integers of the gauge_ids
        gauge_int = [int(item) for item in gauge_ids]
        
        # Load the geometries
        vds = gpd.read_file(sub_catch)
        
        # Get the relevant sub catchments
        vds = vds.astype({"value": int})
        vds = vds[vds.value.isin(gauge_int)]
        
        # Create the data mask with only the gauges of interest
        mask = ds.raster.geometry_mask(vds)
        
        # Data array for soil thickness
        da = ds[ST_key]

        da.values[mask] *= _ST             
        ds[ST_key] = da

        ds.to_netcdf(sml)
        ds.close()
        os.remove(temp_base_grid)
        l.info(f"Writing dataset to {sml}")
    
    #make sure the lakes h-q rating relation is also copied
    lakes_out = [Path(base_dir, 'instates', lakefile.split('/')[-1]) for lakefile in lakes_in]
    if len(lakes_in)!=len(lakes_out):
        l.error(f"lakes_in and lakes_out should have the same length")
        raise ValueError("lakes_in and lakes_out should have the same length")
    for src, dst in zip(lakes_in, lakes_out):
        if not os.path.exists(dst):
            shutil.copy(src, dst)
            l.info(f"Copying {src} to {dst}")

if __name__ == "__main__":
    
    try:
        if "snakemake" in globals():
            mod = globals()["snakemake"]
            l=setup_logging('data/0-log', f"ST_instate_staticmaps_L{mod.params.level}.log")
            main(
                l,
                mod.params.base_grid,
                mod.params.level,
                mod.params.ST_str,
                mod.params.ST_key,
                mod.params.graph[f"level{mod.params.level}"]["elements"],
                mod.params.sub_catch,
                mod.params.lakes_in,
            )

        else:
                # main(
                #     "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/staticmaps.nc",
                #     {"ksat": 5, "rd": 0.5, "st": 0.5},
                #     ["KsatHorFrac", "RootingDepth", "SoilThickness"],
                #     ["set", "mult", "mult"],
                #     ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000'],
                #     "p:/1000365-002-wflow/tmp/usgs_wflow/models/TEST_MODEL_KING/staticgeoms/subcatch_obs.geojson",
                #     "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/calib_data/level1/out.nc"
                # )
            pass

    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e

