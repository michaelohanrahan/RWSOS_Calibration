import os
import sys
import geopandas as gpd
from hydromt.log import setuplog
from hydromt_wflow import WflowModel
from pathlib import Path
from icecream import ic

#python meuse/src/pre/update_gauges.py p:/11209265-grade2023/wflow/RWSOS_Calibration meuse/data/1-external meuse/data/2-interim/QGIS/hourlygauges_snapped.geojson --new_root meuse/data/2-interim/addhourly --basename hourly --index_col id --max_dist 1000 --crs EPSG:4326 --config_fn wflow_sbm_template.toml

os.chdir("p:/11209265-grade2023/wflow/RWSOS_Calibration")
root = os.getcwd()

new_root = "meuse/data/2-interim/addhourly"
basename = "hourly"
index_col = "id"
max_dist = 1000
crs = 'epsg:4326'
derive_subcatch = True
config_fn = "wflow_sbm_template.toml"
gauges = "meuse/data/2-interim/QGIS/hourlygauges_snapped.geojson"
snap_to_river = True
mode = "w+"

if new_root is None:
    new_root = root
    mode = "w+"

logger = setuplog("build", log_level=20)

if not Path(gauges).is_absolute():
    gauges = Path(Path.cwd(), gauges)

w = WflowModel(
    root=os.path.join(root, 'meuse/data/1-external'),
    mode="r",
    config_fn = config_fn,
    data_libs = [],
    logger=logger,
    )

w.read_grid()

w.set_root(
    root=new_root,
    mode=mode
    )

w.set_crs(crs)
ic(w.crs)
# Updating
# Based on gauges geojson
gauges = gpd.read_file(gauges, crs=crs)
if 'geometry' not in gauges.columns:
    if 'x' in gauges.columns and 'y' in gauges.columns:
        gauges['geometry'] = gpd.points_from_xy(gauges['x'], gauges['y'])
    elif 'lon' in gauges.columns and 'lat' in gauges.columns:
        gauges['geometry'] = gpd.points_from_xy(gauges['lon'], gauges['lat'])
    else:
        raise ValueError('Gauges file does not contain geometry, x/y or lon/lat columns')
    
gauges = gauges.set_geometry(col='geometry')
ic(gauges.crs)

w.setup_gauges(
    gauges_fn=gauges,
    snap_to_river=snap_to_river,
    derive_subcatch=derive_subcatch,
    index_col=index_col,
    basename=basename,
    max_dist=max_dist,
    toml_output='csv',
    gauge_toml_header = ["Q", "P", "PET"],
    gauge_toml_param = ["lateral.river.q_av", "vertical.precipitation", "vertical.potential_evaporation"],

)

w.config['input']['path_static'] = 'staticmaps/staticmaps_hourlygauges.nc'
print('writing config')
w.write_config()
print('writing grid')
w.write_grid()
print('writing geoms')
w.write_geoms()