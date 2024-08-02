#!/bin/bash
#SBATCH --job-name=RWSorchestrator                                             # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.log # Standard output and error log
#SBATCH --error=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse/data/0-log/h7/calib_orch_%j.err  # Standard output and error log
#SBATCH --time=0-00:30:00  # Job duration (hh:mm:ss)
#SBATCH --partition test
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Set variables
# Go one directory up to set PWD
cd ..
echo "Changed to $(pwd)"
pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --profile "h7/" --unlock
echo "Unlocked directory"
# echo "Making rulegraph"
# pixi run snakemake -s "2_Snakefile.smk" --configfile "config/calib.yml" --rulegraph | dot -Tsvg > dag/rulegraph.svg
echo "Running snakemake with 1 core"
pixi run snakemake -s "2_Snakefile.smk" -c 1 --configfile "config/calib.yml" --profile "h7/" --forceall --rerun-triggers mtime

