cd ..
echo "Running get_valid_files.py for meuse level5"
pixi run python src/calib/get_valid_files.py "meuse" "level5"
echo "Unlocking meuse level5"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --omit-from $(cat ./valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5  --unlock
echo "Unlocked directory"
echo "cleaning up metadata"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --cleanup-metadata $(cat ./valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5 -n
echo "-n flagged for dry run    *** omitting certain files ***"
pixi run snakemake -s "2_Snakefile.smk" --profile "slurm/" --omit-from $(cat ./valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5 -n
echo "Running snakemake orchestrator only for the remainder of L5 *** omitting certain files ***"
pixi run snakemake -s "2_Snakefile.smk" -c 4 --profile "slurm/" -R wflow_L5 --omit-from $(cat ./valid_outputs.txt | tr ',' ' ') -R wflow_L5  --until done_L5
