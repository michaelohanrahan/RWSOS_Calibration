@echo off

set work_dir=c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine

set cwd=%work_dir%
set new_root=%work_dir%\data\2-interim\gaugesadded
set config_root=%work_dir%\data\1-external
set gauges=%work_dir%\data\2-interim\gaugesadded\all_stations.geojson
set ignore_list=%work_dir%\data\2-interim\ignore_list.txt
set mode=w+
set basename=Hall
set index_col=wflow_id
set max_dist=1000
set crs=EPSG:4326
set config_old=wflow_sbm_template.toml
set config_new=wflow_sbm_addgauges.toml

echo "Running update_gauges.py"
pixi run python c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\src\pre\update_gauges.py "%cwd%" "%config_root%" "%gauges%" --ignore_list "%ignore_list%" --new_root "%new_root%" --mode "%mode%" --basename "%basename%" --index_col "%index_col%" --snap_to_river True --max_dist %max_dist%  --derive_subcatch True --crs "%crs%" --config_old "%config_old%" --config_new "%config_new%"

echo "Running create_dependency_graph.py"
pixi run python c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\src\pre\create_dependency_graph.py "%cwd%" --gridfile "%new_root%\staticmaps\staticmaps.nc" --gaugeset "%basename%" --testmode True

pause