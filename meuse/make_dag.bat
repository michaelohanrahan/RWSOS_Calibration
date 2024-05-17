@echo off

echo Activating the pixi environment
pixi run snakemake -s "snakefile" --configfile "config/calib.yml" --dag | dot -Tsvg > dag.svg

pause