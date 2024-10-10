@echo off

echo Dry Run for the calibration of the RWSOS model final stage
cd ..
pixi run snakemake -s "3_Snakefile_final_stage.smk" --dry-run --configfile "config/calib.yml"

pause