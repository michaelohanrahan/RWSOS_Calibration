@echo off

echo Activating the pixi environment
cd ..
@REM pixi run snakemake -s "2_Snakefile.smk" --unlock --configfile "config/calib.yml" 
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --nolock --rulegraph | dot -Tsvg > dag/rulegraph.svg
pause