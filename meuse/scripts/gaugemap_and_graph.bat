@echo off

set cwd=p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse
set new_root=p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/gaugesadded
set config_root=p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/1-external
set gauges=p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/QGIS/to_wflow/hourly_gauges_all.geojson
set mode=w+
set basename=Hall
set index_col=wflow_id
set max_dist=1000
set crs=EPSG:4326
set config_old=wflow_sbm_template.toml
set config_new=wflow_sbm_addgauges.toml
echo "Running update_gauges.py"
pixi run python src/pre/update_gauges.py %cwd% %config_root% %gauges% --new_root "%new_root%" --mode "%mode%" --basename "%basename%" --index_col "%index_col%" --snap_to_river True --max_dist %max_dist%  --derive_subcatch True --crs "%crs%" --config_old "%config_old%" --config_new "%config_new%"
echo "Running create_dependency_graph.py"
pixi run python src/graph/create_dependency_graph.py %cwd% --gridfile "%new_root%/staticmaps/staticmaps.nc" --gaugeset "%basename%" --testmode True
:: pixi run python src/graph/create_dependency_graph.py p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse --gridfile "p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/gaugesadded/staticmaps/staticmaps.nc" --gaugeset "Hall" --testmode True

pause