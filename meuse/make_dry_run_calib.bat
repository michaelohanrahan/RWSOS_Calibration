@echo off

echo Activating the pixi environment
pixi run snakemake -s "snakefile" --dry-run --configfile "config/calib.yml"

pause