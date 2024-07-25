@echo off

echo Activating the pixi environment
REM Path(DRIVE, '11209265-grade2023', 'wflow', 'RWSOS_Calibration', "meuse"))
cd "p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse"
pixi run snakemake -s "snakefile" --dry-run --configfile "config/calib.yml"

pause