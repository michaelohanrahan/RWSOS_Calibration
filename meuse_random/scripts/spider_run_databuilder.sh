#!/bin/bash
#SBATCH --job-name=data_builder  # Job name
#SBATCH --output=databuilder_%j.log # Standard output and error log log
#SBATCH --time=0-12:00:00  # Job duration (hh:mm:ss)
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p normal
#SBATCH --get-user-env

pixi_toml=/project/afrijnmaas/Data/RWSOS_Calibration/pixi.toml
snakefile=1_Snakefile_databuilder_spider.smk

cd ..
pixi run --manifest-path $pixi_toml snakemake -s $snakefile --unlock
pixi run --manifest-path $pixi_toml snakemake -s $snakefile -n
pixi run --manifest-path $pixi_toml snakemake -s $snakefile