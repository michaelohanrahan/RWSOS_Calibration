from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
from hydromt import DataCatalog


def main():
    pass


if __name__ == "__main__":
    d = DataCatalog(
        data_libs=[
            "deltares_data",
            "p:/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml"
        ]
    )
    geom = gpd.read_file("p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CLIMATE/staticgeoms/region.geojson")
    geom2 = gpd.read_file("p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_PIERCE_CLIMATE/staticgeoms/region.geojson")
    geom = pd.concat([geom, geom2])
    geom = geom.dissolve()
    
    era5 = d.get_rasterdataset("era5_hourly_zarr", time_tuple=("1985-01-01", "2014-12-31T21:00:00"), geom=geom, buffer=2)
    precip = era5.precip
    precip = precip.raster.zonal_stats(geom, "mean")
    cl = d.get_rasterdataset("cmip6_EcEarth_historic", time_tuple=("1985-01-01", "2014-12-31T21:00:00"), geom=geom, buffer=2)
    pass
