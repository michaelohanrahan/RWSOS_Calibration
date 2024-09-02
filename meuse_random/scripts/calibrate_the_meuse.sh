#!/bin/bash
#SBATCH --job-name=RWS_M_Rand_3Kx6L  # Job name
#SBATCH --output=/u/ohanrah/documents/RWSoS/RWSOS_Calibration/meuse_random/data/0-log/h7/RWS_M_Rand_3Kx6L_%j.log # Standard output and error log log
#SBATCH --time=30-00:00:00  # Job duration (hh:mm:ss)
#SBATCH --partition 4pcpu
#SBATCH --ntasks=1  # Number of tasks (analyses) to run
#SBATCH --mail-user=michael.ohanrahan@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

sleep 600

cd ..
echo "Changed to $(pwd)"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --unlock
echo "Unlocked directory"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --dry-run
echo "Running snakemake orchestrator"
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "slurm/" --forceall

