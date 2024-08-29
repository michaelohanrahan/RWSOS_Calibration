#!/bin/bash
cd ..
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock
pixi run snakemake -s "2_Snakefile.smk" -c 1 --configfile "config/calib.yml" -n
pixi run snakemake -s "2_Snakefile.smk" -c 1 --profile "$HOME/.config/snakemake/slurm/"
