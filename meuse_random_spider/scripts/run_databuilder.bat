cd ..
pixi run snakemake -s "1_Snakefile_databuilder.smk" --unlock
pixi run snakemake -s "1_Snakefile_databuilder.smk" -n
pixi run snakemake -s "1_Snakefile_databuilder.smk" --forceall