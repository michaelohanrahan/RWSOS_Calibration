#!/bin/bash
#cd /p/
#Secho "changed to $(pwd)"
cd /p/11209265-grade2023/wflow/RWSOS_Calibration
echo "changed to $(pwd)"
cd meuse
echo "changed to $(pwd)"
echo "Running snakemake with 1 core"
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" -c 1 --nolock --wait-for-files --rerun-incomplete --forceall
