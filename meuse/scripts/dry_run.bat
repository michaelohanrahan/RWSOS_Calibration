@echo off

echo Dry Run for the calibration of the RWSOS model
REM Path(DRIVE, '11209265-grade2023', 'wflow', 'RWSOS_Calibration', "meuse"))
cd ..
pixi run snakemake -s "2_Snakefile.smk" --dry-run --configfile "config/calib.yml" --nolock 

pause