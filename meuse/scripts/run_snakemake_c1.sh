#!/bin/bash
cd ..
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --unlock
HOME = /u/ohanrah/
pixi run snakemake -s "2_Snakefile.smk" --profile $HOME/.config/snakemake/h7/ --configfile "config/calib.yml" -c 1 --nolock --forceall --rerun-incomplete
