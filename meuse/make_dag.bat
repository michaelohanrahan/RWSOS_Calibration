@echo off

echo Activating the pixi environment
pixi run snakemake -s "snakefile" --configfile "config/snakeConfig.yaml" --dag | dot -Tpdf > dag.pdf 

pause