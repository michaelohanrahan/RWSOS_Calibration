# Workflow for the calibration of the wflow model
# This workflow is used to calibrate the wflow model using the data built in the databuilder subworkflow
#@michaelohanrahan
#@jingdeng

import os 
import platform
from pathlib import Path
import yaml

if platform.system() == "Windows":
    DRIVE = "p:/"
    PLATFORM = "Windows"

elif platform.system() == "Linux":
    DRIVE = "/p"
    PLATFORM = "Linux"

configfile: str(Path("config", "calib.yml").as_posix())
#Base directory
basin = config["basin"]    
base_dir = config["base_dir"]  # RWSOS_Calibration
base_dir = f"{DRIVE}/{Path(base_dir).as_posix()}" # p: or /p/ ... / RWSOS_Calibration / basin
workdir: str(Path(base_dir, basin).as_posix())

import json
import shutil
import os
from snakemake.utils import Paramspace
import geopandas as gpd
from src.calib.latin_hyper_paramspace import create_set_all_levels
import glob
import pandas as pd 
import argparse

#=========================================================
log_dir = Path(base_dir, basin, config["log_dir"]).as_posix()          # data/0-log
source_dir = Path(base_dir, basin, config["source_dir"]).as_posix()    # data/1-external
inter_dir = Path(base_dir, basin, config["inter_dir"]).as_posix()      # data/2-interim
input_dir = Path(base_dir, basin, config["input_dir"]).as_posix()      # data/3-input
out_dir = Path(base_dir, basin, config["output_dir"]).as_posix()       # data/4-output
vis_dir = Path(base_dir, basin, config["vis_dir"]).as_posix()          # data/5-visualisation

for dir in [source_dir, inter_dir, out_dir, vis_dir, input_dir]:
    os.makedirs(dir, exist_ok=True)

gauges = config["gauges"]


#=========================================================
#find last level from the final level directory
levels = glob.glob(str(Path(inter_dir,'calib_data', "level*")))
levels_ints = [int(level.split("level")[-1]) for level in levels]
last_level = int(config.get("level", max(levels_ints)))

#define elements from the staticgeoms
elements = list(gpd.read_file(Path(input_dir,"staticgeoms", f'subcatch_{config["gauges"]}.geojson'))["value"].values)

#0:N_samples for each level
N_samples = config["N_SAMPLES"]

#parameter set dataframe using LHS method
lnames, methods, all_level_df = create_set_all_levels(last_level=max(levels_ints), RECIPE=config["calib_recipe"], N_SAMPLES=N_samples, OPTIM='random-cd')
all_level_df.to_csv(Path(inter_dir, "calib_data", "all_level_paramspace.csv"), index=False)

for level in range(0, last_level+1):
    df = all_level_df[all_level_df.index == level]
    df.to_csv(Path(inter_dir, "calib_data", f"level{level}", "paramspace.csv"))

#staticmaps
staticmaps = Path(input_dir, "staticmaps", "staticmaps.nc")

#lakes 
lakes = glob.glob(str(Path(input_dir, "staticmaps", "lake*.csv")))
lakes = [Path(lake).as_posix() for lake in lakes]
lakefiles = [lake.split("/")[-1] for lake in lakes]
assert len(lakes) > 0, "No lakes found, if this is expected comment out this assert"

#subcatch
subcatch = Path(input_dir, "staticgeoms", f'subcatch_{config["gauges"]}.geojson')

#cfg template
cfg_template = Path(input_dir, "wflow_sbm.toml")

#calib_dir
calib_dir = Path(inter_dir, "calib_data")

#graph level
graph = json.load(open(Path(inter_dir, f'{config["gauges"]}_levels_graph.json')))

#graph node
graph_node = json.load(open(Path(inter_dir, f'{config["gauges"]}_nodes_graph.json')))

#graph pred
graph_pred = json.load(open(Path(inter_dir, f'{config["gauges"]}_pred_graph.json')))

#recipe
with open(config["calib_recipe"]) as recipe:
    recipe = json.load(recipe)
snames = [recipe[key]["short_name"] for key in recipe.keys()]

############################
# DOING the snakey!
############################

'''
Preparing the final stage: This rule prepares the final stage of the calibration process.
'''
rule prep_final_stage:
    input: 
        done = Path(calib_dir, "level"+f'{last_level}', "level.done"),
        performance = glob.glob(str(Path(calib_dir, "level*", "performance.zarr"))) #expand(Path(calib_dir, "{level}", "performance.nc"), level=list(graph.keys()))
    params:
        cfg_template = cfg_template,
        cfg_args = [config["eval_runstart"], config["eval_runend"], config["timestep"], Path(source_dir, config["source_forcing_data"])],
        staticmaps = staticmaps
    output: 
        cfg = Path(input_dir, config["wflow_cfg_name"]),
        performance = Path(out_dir, "performance.nc"),
        staticmaps = Path(input_dir, "staticmaps.nc")
    localrule: True
    script:
        """src/calib/prep_final_stage.py"""

'''Final instate: This rule creates the final instate for the model evaluation.'''

rule final_instate_toml:
    input: 
        performance = Path(out_dir, "performance.nc"),
        config_fn = cfg_template,
    params: 
        root = Path(input_dir, "instates").as_posix(),
        level = "final",
        starttime = config["eval_instart"],
        endtime = config["eval_inend"],
        staticmaps = staticmaps
    localrule: True
    output:
        cfg = Path(input_dir, "instates", "post_calib_instate.toml"),
        

rule run_instate:
    input:
        cfg = Path(input_dir, "instates", "post_calib_instate.toml"),
        staticmaps = Path(input_dir, "staticmaps.nc")
    params: 
        project = Path(base_dir, "bin").as_posix(),
    output: 
        outstate=Path(input_dir, "instates", "instate_level_final.nc"),
        done = touch(Path(input_dir, "instates", "done_final_instate.txt"))
    threads: config["wflow_threads"]
    localrule: False
    group: "wflow"
    shell:
        f"""julia --project="{{params.project}}" -t {{threads}} -e \
        "using Pkg;\
        Pkg.instantiate();\
        using Wflow;\
        Wflow.run()" {{input.cfg}}"""

rule run_final_model:
    input:
        done = Path(input_dir, "instates", "done_final_instate.txt"),
        instate = Path(input_dir, "instates", "instate_level_final.nc"),
        cfg = Path(input_dir, config["wflow_cfg_name"]),
        staticmaps = Path(input_dir, "staticmaps.nc")
    params: 
        project = Path(base_dir, "bin").as_posix(),
    output: Path(out_dir, "output_scalar.nc")
    threads: config["wflow_threads"]
    localrule: False
    group: "wflow"
    shell:
        f"""julia --project="{{params.project}}" -t {{threads}} -e \
        "using Pkg;\
        Pkg.instantiate();\
        using Wflow;\
        Wflow.run()" {{input.cfg}}"""

rule visualize:
    input: 
        scalar = Path(out_dir, "output_scalar.nc"),
        performance = Path(out_dir, "performance.nc")
    params:
        observed_data = config["observed_data"],
        gauges = elements,
        starttime = config["eval_starttime"],
        endtime = config["eval_endtime"],
        period_startdate = config["hydro_period_startdate"],
        period_length = config["hydro_period_length"],
        period_unit = config["hydro_period_unit"],
        output_dir = Path(vis_dir, "figures")
    localrule: True
    output:
        figures = expand(Path(vis_dir, "hydro_gauge", "hydro_{gauge}.png"), gauge=elements)
    script:
        """src/post/plot_final_model.py"""