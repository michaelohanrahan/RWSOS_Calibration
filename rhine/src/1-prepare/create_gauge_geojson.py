from pathlib import Path
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point


def gauge_geojson(fn_obs):
    # Create a geojson file for the gauges
    obs = xr.open_dataset(fn_obs)
    
    stations_data = []
    
    for station in obs.wflow_id.values:
        wflow_id = station.item()
        x = obs.sel(wflow_id=station).x.item()
        y = obs.sel(wflow_id=station).y.item()
        
        point = Point(x, y)
        
        stations_data.append({'wflow_id': wflow_id, 'geometry': point})
        
    gdf = gpd.GeoDataFrame(stations_data, crs="EPSG:4326")
    
    gdf.to_file(work_dir / "data/2-interim/gaugesadded/all_stations.geojson", driver='GeoJSON')
        
  
if __name__ == '__main__':
    
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine')
    
    fn_obs = work_dir / 'data/1-external/discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc'
    
    gauge_geojson(fn_obs)