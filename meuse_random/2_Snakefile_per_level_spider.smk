# Workflow for the calibration of the wflow model
# This workflow is used to calibrate the wflow model using the data built in the databuilder subworkflow
#@michaelohanrahan
#@jingdeng

import os
import platform
from pathlib import Path
import yaml

# Spider specific stuff
DRIVE = "/project"
PLATFORM = "Linux"

configfile: str(Path("config", "calib_spider.yml").as_posix())
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
# Container and binding information
wflow_container = "/project/afrijnmaas/Software/containers/julia_container/wflow_julia_v073.sif"

base_toml = "/project/afrijnmaas/Data/RWSOS_Calibration/meuse_random/wflow_sbm_bind.toml"
forcing_fn = "/project/afrijnmaas/Data/RWSOS_Calibration/meuse_random/data/3-input/forcing_Meuse_20050101_20180222_v2_wgs2_remapbil_semisstonn.nc"
#=========================================================

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
:: rule all ::
Defines the overall default target
We expect upon successfull completion that all gauges will have been visualized
'''
rule all:
    input:
        expand(Path(inter_dir, "calib_data", "level{level}", "best_params.csv"), level=range(-1, last_level+1)),
        expand(Path(inter_dir, "calib_data", "level{level}", "performance.zarr"), level=range(0, last_level+1)),
        expand(Path(calib_dir, f"level{level}", "level.done"), level=range(-1, last_level+1))

#THIS RULE ALLOWS RECURSIVE DEPENDENCY
rule init_done:
    params:
        d = Path(calib_dir, "level-1")
    output:
        p = Path(calib_dir, "level-1", "level.done"),
        r = Path(calib_dir, "level-1", "best_params.csv")
    localrule: True
    run:
        os.makedirs(params.d, exist_ok=True)
        with open(output.p, "w") as f:
            f.write("done")
        with open(output.r, "w") as f:
            f.write("")


for _level in range(last_level, last_level+1):
    '''
    ** WILDCARDING **
        -- The wildcarding is done for the level of the calibration
        -- The params are saved in the level directory as paramspace.csv
        -- The paramspace is then wildcarded in the following rules within this loop
    '''
    params_L = Paramspace(all_level_df[all_level_df.index == _level])
    '''
    :: Random dataframe ::
        SUMMARY??
    '''
    rule:
        name: f"random_data_L{_level}"
        input:
            done = Path(calib_dir, f"level{_level-1}", "level.done"),
            best_params_previous = Path(calib_dir, f"level{_level-1}", "best_params.csv")
        params:
            level = f"level{_level}",
            params_df = df,
            graph = graph,                      # Hall_levels_graph.json
            graph_pred = graph_pred,            # Hall_pred_graph.json
            graph_node = graph_node,            # Hall_nodes_graph.json
        # localrule: True
        output:
            random_params = Path(calib_dir, f"level{_level}", "random_params.csv")
        script:
            """src/calib/random_params.py"""

    """
    ::config::
            This rule modifies a blueprint configuration file for a specific time period and forcing path.
            It takes the blueprint configuration file, start time, end time, time step, and forcing path as parameters.
            The output is a set of modified configuration files.
    """
    rule:
        name: f"config_L{_level}"
        input:
            random_params = Path(calib_dir, f"level{_level}","random_params.csv")
        params:
            level = _level,
            cfg_template = Path(cfg_template).as_posix(),
            starttime = config["starttime"],
            endtime = config["endtime"],
            forcing_path = Path(input_dir, config["source_forcing_data"]),
            gaugemap=f"gauges_{config['gauges']}",
            wildness = params_L.wildcard_pattern,
        # localrule: True
        output:
            out_file=Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, config["wflow_cfg_name"])
        script:
            """src/calib/set_config.py"""

    '''
    :: create_params ::
            This rule creates parameter sets for the Wflow model.
            It takes a configuration file, a dataset of static maps, a parameter space instance,
            a list of parameter names, a list of parameter methods, a level, a graph, and a sub-catchment as parameters.
            The output is a set of static map files.
    NEW: The parameter space is wildcarded but if upstream bestparams exist we will randomly choose which one to apply.
    '''

    rule:
        name: f"create_params_L{_level}"
        input:
            config_fn = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, config["wflow_cfg_name"]), #wait for the config to be done
            random_params = Path(calib_dir, f"level{_level}", "random_params.csv"),
        params:
            dataset = staticmaps,
            params = params_L.instance,
            params_sname = snames,
            params_lname = lnames,
            params_method = methods,
            level = f"level{_level}",
            graph = graph,
            sub_catch = subcatch,
            lake_in = lakes
        # localrule: True
        output:
            staticmaps = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "staticmaps.nc"),
            # model_dir = directory(Path(calib_dir, f"level{_level}", params_L.wildcard_pattern)),
            lake_out = expand(Path(calib_dir, f"level{_level}", "{params}", "{lakes}"), params=params_L.wildcard_pattern, lakes=lakefiles),
        script:
            """src/calib/set_calib_params.py"""

    '''
    ::wflow::
            This rule runs the Wflow model. It takes a configuration file and a grid (staticmap) file as input.
    '''
    rule:
        name: f"wflow_L{_level}"
        input:
            cfg = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, config["wflow_cfg_name"]),
            staticmaps = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "staticmaps.nc"),
            lake_hqs = expand(Path(calib_dir, f"level{_level}", "{params}", "{lakes}"), params=params_L.wildcard_pattern, lakes=lakefiles),
            # model_dir = directory(Path(calib_dir, f"level{_level}", params_L.wildcard_pattern))
        params:
            project = Path(base_dir, "bin").as_posix(),
        output:
            Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "output_scalar.nc"),
        # localrule: False
        # group: f"wflow_L{_level}"
        resources:
            time = "12:00:00",
            mem_mb = 8000
        run:
            shell(
                """
                apptainer run \
                    --bind {tmp_model_dir}:/data \
                    --bind {tmp_base_toml}:/data/wflow_sbm.toml \
                    --bind {tmp_forcing_fn}:/data/inmaps.nc \
                    --writable-tmpfs \
                    {wflow_container} "using Wflow; Wflow.run()" "/data/wflow_sbm.toml"
                """.format(tmp_base_toml = base_toml, tmp_forcing_fn=forcing_fn, tmp_model_dir=Path(input.staticmaps).parent, wflow_container=wflow_container)
            )
    rule:
        name: f"done_L{_level}"
        input:
            output = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "output_scalar.nc")
        output:
            done = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "run.done")
        localrule: True
        shell:
            """touch {output.done}"""

    '''
    Evaluate: This rule evaluates the performance of the model by comparing the model output with observed data and outputs the best parameters and performance metrics.
    Output:  an unstacked performance.nc file with metrics and weights for each parameter set, and a best_10params.csv file with the best parameter set.
            out csv with the best parameter set for each gauge

    '''
    rule:
        name: f"evaluate_L{_level}"
        input:
            sim = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "output_scalar.nc"),
            done = expand(Path(calib_dir, f"level{_level}", "{params}", "run.done"), params=params_L.instance_patterns),
        params:
            observed = Path(source_dir, config["observed_data"]),
            dry_month = config["dry_months"],
            window = config["window"],
            level = f"{_level}",
            graph = graph,
            params = params_L.wildcard_pattern,
            starttime = config["cal_eval_starttime"],
            endtime = config["cal_eval_endtime"],
            metrics = config["metrics"], #["kge", "nselog_mm7q", "mae_peak_timing", "mape_peak_magnitude"]
            weights = config["weights"],
            gaugeset = f"Q_gauges_{config['gauges']}",
        # localrule: False
        # group: f"evaluate_L{_level}"   #TODO: add this group to the config
        output:
            performance = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "performance.nc"),
            eval_done = Path(calib_dir, f"level{_level}", params_L.wildcard_pattern, "evaluate.done"),
        resources:
            time = "00:30:00",
            mem_mb = 8000
        script:
            """src/calib/evaluate_per_run.py""" # TODO: modify this script to submit multiple grouped jobs

    rule:
        name: f"combine_performance_L{_level}"
        input:
            done = expand(Path(calib_dir, f"level{_level}", '{params}', "evaluate.done"), params=params_L.instance_patterns),
            performance_files=expand(Path(calib_dir, f"level{_level}", "{params}", "performance.nc"), params=params_L.instance_patterns)
        output:
            performance = Path(calib_dir, f"level{_level}", "performance.zarr"),
            best_params = Path(calib_dir, f"level{_level}", "best_params.csv") #defaults to best 10
        resources:
            time = "01:00:00",
            mem_mb = 32000
        script:
            "src/calib/combine_evaluated.py"

    '''
    :: set params ::
            This rule inherits the best parameters from the previous level and sets those within the base staticmaps.
            The previous level staticmaps is saved for posterity in the intermediate folder.
            The output is a set of staticmaps and a level.done file.
    '''

    rule:
        name: f"set_params_L{_level}"
        input:
            best_params = Path(calib_dir, f"level{_level}", "best_params.csv"),
        params:
            staticmaps = staticmaps,
            sub_catch = subcatch,
            params_lname = lnames,
            params_method = methods,
            level=_level
        # localrule: True
        output:
            done_nc = Path(input_dir, "staticmaps", "intermediate", f"staticmaps_L{_level-1}.nc"),
            done = Path(calib_dir, f"level{_level}", "level.done")
        script:
            """src/calib/set_eval_params.py"""

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
    # localrule: True
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
    # localrule: True
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
    # localrule: True
    output:
        figures = expand(Path(vis_dir, "hydro_gauge", "hydro_{gauge}.png"), gauge=elements)
    script:
        """src/post/plot_final_model.py"""