from pathlib import Path
import xarray as xr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Point


def main(root: str,
         observations: str,
         selected_gauge: str,
         outdir: str):
    # Create a geojson file for the gauges
    work_dir = Path(root)
    fn_obs = work_dir / observations
    fn_selected_gauge = work_dir / selected_gauge
    
    obs = xr.open_dataset(fn_obs)
    df = pd.read_csv(fn_selected_gauge)
    
    # extract numbers (gauge wflow_id from dataframe)
    _temp = df.values.flatten()
    sel_gauges = [int(num) for num in _temp if pd.notna(num)]
    
    stations_data = []
    
    for station in sel_gauges:
        wflow_id = station
        if station == 709:
            # the original location of Lobith (709) is outside of model boundary
            # thus, manually assigned it to be inside of boundary
            x = 6.106169988948045
            y = 51.849728794940447
        else:
            x = obs.sel(wflow_id=station).x.item()
            y = obs.sel(wflow_id=station).y.item()
        
        point = Point(x, y)
        
        stations_data.append({'wflow_id': wflow_id, 'geometry': point})
        
    gdf = gpd.GeoDataFrame(stations_data, crs="EPSG:4326")
    
    gdf.to_file(work_dir / outdir, driver='GeoJSON')
        
  
if __name__ == '__main__':
    
    root = r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine'
    observations = r'data/1-external/discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc'
    selected_gauge = r'data/2-interim/manually_selected_gauges.csv'
    outdir = r"data/2-interim/gaugesadded/all_stations.geojson"
    
    main(root,
         observations,
         selected_gauge,
         outdir)
    
