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

os.chdir(Path(DRIVE, "11209265-grade2023", "wflow", "RWSOS_Calibration", "meuse_random_spider"))

configfile: str(Path("config", "calib.yml").as_posix())
#Base directory
basin = config["basin"]    
base_dir = config["base_dir"]  # RWSOS_Calibration
base_dir = f"{DRIVE}/{Path(base_dir).as_posix()}" # p: or /p/ ... / RWSOS_Calibration / basin
workdir: str(Path(base_dir, basin).as_posix())

import json
import shutil
import os
import snakemake
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
# all_level_df.to_csv(Path(inter_dir, "calib_data", "all_level_paramspace.csv"), index=False)

# for level in range(0, last_level+1):
#     df = all_level_df[all_level_df.index == level]
#     df.to_csv(Path(inter_dir, "calib_data", f"level{level}", "paramspace.csv"))

#staticmaps
staticmaps = Path(input_dir, "staticmaps", "staticmaps.nc")

#lakes 
lakes = glob.glob(str(Path(input_dir, "staticmaps", "lake*.csv")))
lakes = [Path(lake).as_posix() for lake in lakes]
lakefiles = [lake.split("/")[-1] for lake in lakes]
# assert len(lakes) > 0, "No lakes found, if this is expected comment out this assert"

#subcatch
subcatch = Path(input_dir, "staticgeoms", f'subcatch_{config["gauges"]}.geojson')

#cfg template
cfg_template = Path(input_dir, "wflow_sbm.toml")

#calib_dir
calib_dir = Path(inter_dir, "calib_data")

# #graph level
# graph = json.load(open(Path(inter_dir, f'{config["gauges"]}_levels_graph.json')))

# #graph node
# graph_node = json.load(open(Path(inter_dir, f'{config["gauges"]}_nodes_graph.json')))

# #graph pred
# graph_pred = json.load(open(Path(inter_dir, f'{config["gauges"]}_pred_graph.json')))

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

# # temporary
# yaml_path = Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random_spider\config\calib.yml')
# with open(yaml_path, 'r') as file:
#     config = yaml.safe_load(file)
    

best_params_fn = Path(inter_dir, "calib_data", "level5", "best_params.csv")
best_params = pd.read_csv(best_params_fn, index_col=['level', 'gauge'])
TOP_ENS = best_params.columns.values

rule all:
    input:
        expand(Path(out_dir, "output_{Topx}", "output_scalar.nc"), Topx=TOP_ENS)


"""Prepare final stage: set staticmaps using the best params"""
rule final_staticmaps: # get parameter sets from Top_x, set parameter values to the staticmaps, prepare toml file
    input:
        best_params_fn = Path(inter_dir, "calib_data", "level5", "best_params.csv")
    params:
        dataset = staticmaps,
        best_params = best_params,
        topx = lambda wildcards: wildcards.Topx,
        params_sname = snames,
        params_lname = lnames,
        params_method = methods,
        sub_catch = subcatch,
        lake_in = lakes
    output: 
        staticmaps = Path(input_dir, "input_{Topx}", "staticmaps.nc"),
        lake_out = expand(Path(input_dir, "input_{Topx}", "{lakes}"), Topx="{Topx}", lakes=lakefiles)
    localrule: True
    script:
        """src/calib/set_staticmaps_final.py"""


'''Final instate: This rule creates the final instate for the model evaluation.'''

rule final_instate_toml:
    input: 
        config_fn = cfg_template,
    params:
        root = Path(input_dir, "input_{Topx}", "instates").as_posix(), 
        starttime = config["eval_instart"],
        endtime = config["eval_inend"],
        forcing_path = Path(input_dir, config["source_forcing_data"]),
        topx = lambda wildcards: wildcards.Topx
    localrule: True
    output:
        cfg = Path(input_dir, "input_{Topx}", "instates", "post_calib_instate.toml")
    script:
        """src/calib/create_instate_final.py"""
        

rule run_instate:
    input:
        cfg = Path(input_dir, "input_{Topx}", "instates", "post_calib_instate.toml"),
        staticmaps = Path(input_dir, "input_{Topx}", "staticmaps.nc")
    params: 
        project = Path(base_dir, "bin").as_posix(),
        topx = lambda wildcards: wildcards.Topx
    output: 
        outstate=Path(input_dir, "input_{Topx}","instates", "instate_final.nc"),
        done = touch(Path(input_dir, "input_{Topx}","instates", "done_final_instate.txt"))
    localrule: False
    threads: 4
    resources:
        mem_mb = 32000,
        time = "5:00:00"
    shell:
        f"""julia --project="{{params.project}}" -t {{threads}} -e \
        "using Pkg;\
        Pkg.instantiate();\
        using Wflow;\
        Wflow.run()" {{input.cfg}}"""


rule final_config_toml:
    input:
        config_fn = cfg_template,
        instate = Path(input_dir, "input_{Topx}","instates", "instate_final.nc")
    params:
        root = Path(input_dir, "input_{Topx}").as_posix(), 
        starttime = config["eval_runstart"],
        endtime = config["eval_runend"],
        forcing_path = Path(input_dir, config["source_forcing_data"]),
        topx = lambda wildcards: wildcards.Topx,
        gaugemap = f"gauges_{config['gauges']}"
    localrule: True
    output:
        cfg = Path(input_dir, "input_{Topx}", config["wflow_cfg_name"])
    script:
        """src/calib/create_config_final.py"""


rule run_final_model:
    input:
        done = Path(input_dir, "input_{Topx}","instates", "done_final_instate.txt"),
        cfg = Path(input_dir, "input_{Topx}", config["wflow_cfg_name"]),
        staticmaps = Path(input_dir, "input_{Topx}","staticmaps.nc")
    params: 
        project = Path(base_dir, "bin").as_posix(),
        topx = lambda wildcards: wildcards.Topx
    output: Path(out_dir, "output_{Topx}", "output_scalar.nc")
    localrule: False
    threads: 4
    resources:
        mem_mb = 32000,
        time = "10:00:00"
    shell:
        f"""julia --project="{{params.project}}" -t {{threads}} -e \
        "using Pkg;\
        Pkg.instantiate();\
        using Wflow;\
        Wflow.run()" {{input.cfg}}"""


# TODO: to be modified: add eval metrics+signatures
# rule final_eval_metrics:
# """dataset of all eval metrics for final models: kge, nse, nse_log, nse_mm7q, mae_pt, mape_pm
# dim: runs (obs, hbv, Topx), wflow_id
# """


rule visualize_one:
    input: 
        scalar = Path(out_dir, "output_{Topx}", "output_scalar.nc")
        # performance = Path(out_dir, "performance.nc")
    params:
        observed_data = config["observed_data"],
        gauges = elements,
        starttime = config["eval_starttime"],
        endtime = config["eval_endtime"],
        topx = lambda wildcards: wildcards.Topx,
        run_list = ["{Topx}"], #This is relatively new 
        period_startdate = config["hydro_period_startdate"],
        period_length = config["hydro_period_length"],
        period_unit = config["hydro_period_unit"],
        output_dir = Path(vis_dir, "figures", "Topx_{Topx}")
    localrule: True
    output:
        figures = Path(vis_dir, "per_run", "Topx_{Topx}", f"hydro_{elements[0]}.png") #cant expand just one wildcard but will create all gauge plots
    script:
        """src/post/plot_final_model.py"""


rule visualize_all:
    input:
        scalar = expand(Path(out_dir, "output_{Topx}", "output_scalar.nc"), Topx=TOP_ENS)
        # performance = Path(out_dir, "performance.nc")
    params:
        observed_data = config["observed_data"],
        gauges = elements,
        starttime = config["eval_starttime"],
        endtime = config["eval_endtime"],
        period_startdate = config["hydro_period_startdate"],
        period_length = config["hydro_period_length"],
        period_unit = config["hydro_period_unit"],
        output_dir = Path(vis_dir, "all_runs")
    localrule: True
    output:
        figures = expand(Path(vis_dir, "all_runs", "hydro_{gauge}.png"), gauge=elements)
        
