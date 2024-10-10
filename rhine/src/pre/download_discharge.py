import matplotlib.pyplot as plt

import logging
import datetime
import numpy as np
import xarray as xr
from tqdm import tqdm
import geopandas as gpd
import hydrofunctions as hf

# Prepare logger
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def add_metadata(meta_dict, site, contrib_area):
    meta = {
        "site_name": meta_dict[f"USGS:{site}"]["siteName"],
        "latitude": meta_dict[f"USGS:{site}"]["siteLatLongSrs"]["latitude"],
        "longitude": meta_dict[f"USGS:{site}"]["siteLatLongSrs"]["longitude"],
        "drainage_area": contrib_area,
        "drainage_area_unit": "square miles",
        "unit": meta_dict[f"USGS:{site}"]["timeSeries"][variable]["variableUnit"],
        "description": meta_dict[f"USGS:{site}"]["timeSeries"][variable][
            "variableDescription"
        ],
    }
    return meta


# Constants
start = "1984-01-01"
end = "2023-01-01"
type = "iv"  # instantaneous values (or dv for daily values)
variable = "00060"  # for discharge, cubic feet per second
resample_freq = "1H"
output_dir = R"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES"
output_file = f"{output_dir}/discharge_obs_hourly.nc"

# Get sites from geopackages
gdf_king = gpd.read_file(
    R"p:\1000365-002-wflow\tmp\usgs_wflow\models\GIS\gauges_king.gpkg"
)
gdf_pierce = gpd.read_file(
    R"p:\1000365-002-wflow\tmp\usgs_wflow\models\GIS\gauges_pierce.gpkg"
)

# Get list of sites
sites = list(gdf_king["SiteNumber"]) + list(gdf_pierce["SiteNumber"])
sites = list(map(str, sites))

# create empty dataset
try:
    ds = xr.open_dataset(output_file).load()
    ds.close()
except:
    ds = xr.Dataset()

# Loop through all sites
for site in tqdm(sites):
    if site not in ds.data_vars:
        # Get information per site
        info_site = hf.site_file(site, verbose=False)
        name_loc = info_site.table.station_nm[0]
        latitude_loc = info_site.table.dec_lat_va[0]
        longitude_loc = info_site.table.dec_long_va[0]
        drain_area = info_site.table.drain_area_va[0]

        # Request data
        request = hf.NWIS(site, type, start, end, verbose=False)
        try:
            logging.info(f"{site} - requesting data")
            df_q = request.df(variable)  # discharge

            # fix index
            index = df_q.index
            index_array = index.to_numpy(dtype=np.datetime64)
            index_array = index_array.astype(datetime.datetime)
            df_q.index = index_array
            df_q.index.name = "time"

            # Resample to 3 hourly data
            df_q_resampled = df_q.resample(resample_freq).mean()

            # To series and correct name
            sr_tmp = df_q_resampled[df_q_resampled.columns[0]]
            sr_tmp.name = site

            # To dataarray
            da_tmp = xr.DataArray.from_series(sr_tmp)
            da_tmp.attrs = add_metadata(
                meta_dict=request.meta, site=site, contrib_area=drain_area
            )
            ds = xr.merge([ds, da_tmp])

        except ValueError:
            logging.warning(f"{site} - variable not found, skipping")
            pass

# Write to file
logging.info(f"{site} - start writing dataset")
comp = dict(zlib=True, complevel=1)
encoding = {var: comp for var in ds.data_vars}
ds.to_netcdf(output_file, encoding=encoding)
