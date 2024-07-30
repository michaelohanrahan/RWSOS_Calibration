@echo off

echo Activating the pixi environment
cd "p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse"
pixi run snakemake -s "Snakefile" --configfile "config/calib.yml" --rulegraph | dot -Tsvg > rulegraph.svg

pause