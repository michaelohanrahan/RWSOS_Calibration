#!/bin/bash
cd ..
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --rulegraph | dot -Tsvg > dag/rulegraph.svg
read -p "Press any key to continue..."