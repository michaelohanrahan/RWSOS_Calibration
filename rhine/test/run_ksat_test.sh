#!/bin/bash
#SBATCH --job-name=ksat_test
#SBATCH --cpus-per-task=1
#SBATCH --partition test
#SBATCH --ntasks=1
#SBATCH --mail-type=all
#SBATCH --time=00:30:00
#SBATCH --mail-user=michael.ohanrahan@deltares.nl


echo "current working directory: $PWD"
echo "calculating PET"

pixi run snakemake -s "/p/11209265-grade2023/wflow/RWSOS_Calibration/rhine/test/ksat_test.smk" --cores 4 --forceall --unlock
pixi run snakemake -s "/p/11209265-grade2023/wflow/RWSOS_Calibration/rhine/test/ksat_test.smk" --cores 4 --forceall -n --quiet rules
pixi run snakemake -s "/p/11209265-grade2023/wflow/RWSOS_Calibration/rhine/test/ksat_test.smk" --cores 4 --profile './test_run_profile.yml'