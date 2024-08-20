@echo off

echo Dry Run for the calibration of the RWSOS model
cd ..
pixi run snakemake -s "2_Snakefile.smk" --dry-run --configfile "config/calib.yml" --nolock --profile "/u/ohanrah/.config/snakemake/slurm/"

pause