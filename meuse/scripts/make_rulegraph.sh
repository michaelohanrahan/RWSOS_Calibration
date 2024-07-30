#!/bin/bash
pixi run snakemake -s "Snakefile" --configfile "config/calib.yml" --rulegraph | dot -Tsvg > dag/rulegraph.svg
read -p "Press any key to continue..."