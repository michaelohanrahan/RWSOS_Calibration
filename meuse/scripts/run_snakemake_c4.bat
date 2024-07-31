@echo off

echo 4 core calibration of the RWSOS model
REM Path(DRIVE, '11209265-grade2023', 'wflow', 'RWSOS_Calibration', "meuse"))
REM cd /d "p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse"
REM echo current directory: %cd%
REM pixi run snakemake -s "2_Snakefile" --configfile "config/calib.yml" unlock
cd ..

pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock
snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" -c 4 --nolock --wait-for-files --rerun-incomplete --forceall
REM scripts/snakemake_c4.bat
pause
