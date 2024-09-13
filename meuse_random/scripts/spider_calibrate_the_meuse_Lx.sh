#!/bin/bash
#SBATCH --time=5-00:00:00  # Job duration (hh:mm:ss)
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --get-user-env

pixi_toml=/project/afrijnmaas/Data/RWSOS_Calibration/pixi.toml
snakefile=2_Snakefile_per_level_spider.smk
profile=/project/afrijnmaas/Data/RWSOS_Calibration/meuse_random/spider

# Go up one directory to ensure correct paths
cd ..

#Unlocking the directory for snakemake
pixi run --manifest-path $pixi_toml snakemake -s $snakefile --unlock --profile $profile --config level=$1
echo "Unlocked directory"
echo "Executing Dry Run"
pixi run --manifest-path $pixi_toml snakemake -s $snakefile --dry-run --quiet --profile $profile --config level=$1

echo "Running snakemake orchestrator"
pixi run --manifest-path $pixi_toml snakemake -s $snakefile --profile $profile --config level=$1

