#!/bin/bash
#SBATCH --job-name=knmi23_Ld          # Job name
#SBATCH --output=/u/buitink/ribasim_rijn/logs/Ld_snakemake_run_%j.log     # Standard output and error log
#SBATCH --time=10-00:00:00           # Job duration (hh:mm:ss)
#SBATCH --partition 1vcpu
#SBATCH --ntasks=1                  # Number of tasks (analyses) to run
#SBATCH --mail-user=joost.buitink@deltares.nl
#SBATCH --mail-type=ALL
#SBATCH --get-user-env

# Initiating snakemake and running workflow in cluster mode
source /u/buitink/miniconda3/bin/activate grade_climate
# conda config --set channel_priority strict

#Going to the folder where scripts are
ROOT="/u/buitink/ribasim_rijn"
cd "${ROOT}"

export name=Ld
export rule=all
export scen="L"

#Unlocking the directory for snakemake
snakemake --unlock -s Snakefile-$scen-v3_newDates --configfile snake_config/snake_settings_$name.yml

#All the cluster configuration (both general and rule specific) is in ~/.config/snakemake/simple
# snakemake $rule -n -q -s Snakefile-$scen --configfile snake_config/snake_settings_Hd.yml --profile knmi23_240year/ --wait-for-files --directory $PWD  --rerun-triggers mtime #--group-components preprocess=3120 xr_merge=50  #--retries 2 --allowed-rules run_wflow
snakemake $rule -s Snakefile-$scen-v3_newDates --configfile snake_config/snake_settings_$name.yml --profile knmi23_240year/ --wait-for-files --directory $PWD  --rerun-triggers mtime #--group-components preprocess=3120 xr_merge=50  #--retries 2 --allowed-rules run_wflow
#snakemake -s snakefile --configfile config/members_config.yml --cluster "sbatch --time=0-00:20:00 --partition 4vcpu --cpus-per-task=4" --jobs 20 --latency-wait 60 --wait-for-files --use-conda --directory $PWD --keep-going --group-components preprocess=3120 xr_merge=10 --allowed-rules run_wflow #--retries 2 --profile simple/

conda deactivate
