@echo off

set work_dir=c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine

set cwd=%work_dir%

set observations=data\1-external\discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc
set selected_gauge=data\2-interim\manually_selected_gauges.csv
set outdir=data\2-interim\gaugesadded\all_stations.geojson

set new_root=%work_dir%\data\2-interim\gaugesadded
set config_root=%work_dir%\data\1-external
set gauges=%work_dir%\data\2-interim\gaugesadded\all_stations.geojson
set ignore_list=None
set mode=w+
set basename=Hall
set index_col=wflow_id
set max_dist=1000
set crs=EPSG:4326
set config_old=wflow_sbm_template.toml
set config_new=wflow_sbm_addgauges.toml

echo "Running create_gauge_geojson.py"
pixi run python c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\src\pre\create_gauge_geojson.py --root "%cwd%" --observations "%observations%" --selected_gauge "%selected_gauge%" --outdir "%outdir%" 

echo "Running update_gauges.py"
pixi run python c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\src\pre\update_gauges.py "%cwd%" "%config_root%" "%gauges%" --ignore_list "%ignore_list%" --new_root "%new_root%" --mode "%mode%" --basename "%basename%" --index_col "%index_col%" --snap_to_river True --max_dist %max_dist%  --derive_subcatch True --crs "%crs%" --config_old "%config_old%" --config_new "%config_new%"

echo "Running create_dependency_graph.py"
pixi run python c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\src\pre\create_dependency_graph.py "%cwd%" --gridfile "%new_root%\staticmaps\staticmaps.nc" --gaugeset "%basename%" --testmode False

pause