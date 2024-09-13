
log_dir=/project/afrijnmaas/Data/RWSOS_Calibration/meuse_random/data/0-log/spider

# Start level0
level0id=$(sbatch --job-name=meuse_level0 --output=$log_dir/snakemake_level0_%j.log spider_calibrate_the_meuse_Lx.sh 0 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level0 with id: $level0id"

# Start level1
level1id=$(sbatch --job-name=meuse_level1 --dependency=afterok:$level0id --output=$log_dir/snakemake_level1_%j.log spider_calibrate_the_meuse_Lx.sh 1 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level1 with id: $level1id"

# Start level1
level2id=$(sbatch --job-name=meuse_level1 --dependency=afterok:$level1id --output=$log_dir/snakemake_level2_%j.log spider_calibrate_the_meuse_Lx.sh 2 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level2 with id: $level2id"

# Start level1
level3id=$(sbatch --job-name=meuse_level1 --dependency=afterok:$level2id --output=$log_dir/snakemake_level3_%j.log spider_calibrate_the_meuse_Lx.sh 3 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level3 with id: $level3id"

# Start level1
level4id=$(sbatch --job-name=meuse_level1 --dependency=afterok:$level3id --output=$log_dir/snakemake_level4_%j.log spider_calibrate_the_meuse_Lx.sh 4 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level4 with id: $level4id"

# Start level1
level5id=$(sbatch --job-name=meuse_level1 --dependency=afterok:$level4id --output=$log_dir/snakemake_level5_%j.log spider_calibrate_the_meuse_Lx.sh 5 | awk 'match($0, /[0-9]+/) {print substr($0, RSTART, RLENGTH)}')
echo "started level5 with id: $level5id"
