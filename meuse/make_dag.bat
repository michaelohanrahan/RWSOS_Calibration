@echo off

echo Activating the pixi environment
<<<<<<< HEAD
pixi run snakemake -s "snakefile" --configfile "config/snakeConfig.yaml" --dag | dot -Tpdf > dag.pdf 
=======
pixi run snakemake -s "snakefile" --configfile "config/calib.yml" --dag | dot -Tsvg > dag.svg
>>>>>>> 98b8f6a20bb35496c5a9fdac156eb997b2436951

pause