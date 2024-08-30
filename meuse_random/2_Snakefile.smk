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
basin = config["basin"]        # meuse
base_dir = config["base_dir"]  # RWSOS_Calibration
base_dir = f"{DRIVE}/{Path(base_dir).as_posix()}" # p: or /p/ ... / RWSOS_Calibration / basin
workdir: str(Path(base_dir, basin).as_posix())
# profile: str(Path(base_dir, basin, "config", "slurm").as_posix())


import json
import shutil
import os
from snakemake.utils import Paramspace
import geopandas as gpd
from src.calib.latin_hyper_paramspace import create_set_all_levels
import glob

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


#find last level from the final level directory
levels = glob.glob(str(Path(inter_dir,'calib_data', "level*")))
levels_ints = [int(level.split("level")[-1]) for level in levels]
last_level = 5 #int(levels[-1].split("level")[-1])

#define elements from the staticgeoms
elements = list(gpd.read_file(Path(input_dir,"staticgeoms", f'subcatch_{config["gauges"]}.geojson'))["value"].values)

#parameter set dataframe using LHS method
#TODO: N_SAMPLES to config
lnames, methods, all_level_df = create_set_all_levels(last_level=last_level, RECIPE=config["calib_recipe"], N_SAMPLES=config["N_SAMPLES"], OPTIM='random-cd')

#0:N_samples for each level
N_samples = config["N_SAMPLES"]


#staticmaps
staticmaps = Path(input_dir, "staticmaps", "staticmaps.nc")

#lakes 
lakes = glob.glob(str(Path(input_dir, "staticmaps", "lake*.csv")))
lakes = [Path(lake).as_posix() for lake in lakes]
lakefiles = [lake.split("/")[-1] for lake in lakes]

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

#soilthicknesses
with open(config["calib_recipe"]) as recipe:
    recipe = json.load(recipe)
    # print(recipe)

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
        expand(Path(vis_dir, "hydro_gauge", "hydro_{gauge}.png"), gauge=elements),
        expand(Path(calib_dir, "level{nlevel}", "level.done"), nlevel=range(-1, last_level+1)),
    default_target: True

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
            f.write()

for _level in range(0, last_level+1):
    
    '''
    ** WILDCARDING **
    -- The wildcarding is done for the level of the calibration
    #¿¿Q??: Does this actually work ot build the dag or do we have to slice the wildcards instead??? 
    '''
    # slice the parameter dataframe for the current level
    #slice 0*2999:1*2999, 3000:5999, 6000:8999, etc
    df = all_level_df.iloc[_level*N_samples-1:(1+_level)*N_samples-1]
    paramspace = Paramspace(df)
    
    '''
    :: Random dataframe ::
        -- In the previous 
    '''
    rule:
        name: f"random_data_L{_level}"
        input: 
            done = Path(calib_dir, f"level{_level-1}", "level.done"),
            best_params_previous = Path(calib_dir, f"level{_level-1}", "best_params.csv")
        params:
            level = f"level{_level}",
            params_df = df,
            graph = graph,  # Hall_levels_graph.json
            graph_pred = graph_pred,  # Hall_pred_graph.json
            graph_node = graph_node,  # Hall_nodes_graph.json
        localrule: True
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
            # flag = Path(calib_dir, f"level{_level-1}","done.txt"),
            random_params = Path(calib_dir, f"level{_level-1}","random_params.csv")
        params: 
            level = _level,
            cfg_template = Path(cfg_template).as_posix(),
            starttime = config["starttime"],
            endtime = config["endtime"],
            forcing_path = Path(input_dir, config["source_forcing_data"]),
            gaugemap=f"gauges_{config['gauges']}",
        localrule: True
        output: 
            expand(Path(calib_dir, f"level{_level}", "{params}", config["wflow_cfg_name"]), params=paramspace.wildcard_pattern)
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
            Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, config["wflow_cfg_name"]),  #TODO: what is this for?
            random_params = Path(calib_dir, f"level{_level}", "random_params.csv")
        params:
            dataset = staticmaps,
            params = paramspace.instance,
            params_lname = lnames,
            params_method = methods,
            level = f"level{_level}",
            graph = graph,
            sub_catch = subcatch,
            lake_in = lakes
        localrule: True
        output: 
            staticmaps = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "staticmaps.nc"),
            lake_out = expand(Path(calib_dir, f"level{_level}", "{params}", "{lakes}"), params=paramspace.wildcard_pattern, lakes=lakefiles)
        script: 
            """src/calib/set_calib_params.py"""

    '''
    ::wflow:: 
            This rule runs the hydrological model. It takes a configuration file and a grid (staticmap) file as input.
    '''
    rule:
        name: f"wflow_L{_level}"
        input: 
            cfg = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, config["wflow_cfg_name"]),
            staticmaps = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "staticmaps.nc"),
            lake_hqs = expand(Path(calib_dir, f"level{_level}", "{params}", "{lakes}"), params=paramspace.wildcard_pattern, lakes=lakefiles)
        params:
            project = Path(base_dir, "bin").as_posix()
        output: 
            Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "output_scalar.nc")
        localrule: False
        group: f"wflow_L{_level}"
        threads: 1
        resources: 
            time = "12:00:00",
            mem_mb = 8000
        shell: 
            f"""julia --project="{{params.project}}" -t {{threads}} -e \
            "using Pkg;\
            Pkg.instantiate();\
            using Wflow;\
             Wflow.run()" {{input.cfg}}"""
    
    rule:
        name: f"done_L{_level}"
        input: 
            output = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "output_scalar.nc")
        output:
            done = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "run.done")
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
            sim = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "output_scalar.nc"),
            done = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "run.done"),
        params: 
            observed = Path(source_dir, config["observed_data"]),
            dry_month = config["dry_months"],
            window = config["window"],
            level = f"{_level}",
            graph = graph,
            params = paramspace.wildcard_pattern, 
            starttime = config["cal_eval_starttime"],
            endtime = config["cal_eval_endtime"],
            metrics = config["metrics"], #["kge", "nselog_mm7q", "mae_peak_timing", "mape_peak_magnitude"]
            weights = config["weights"], 
            gaugeset = f"Q_gauges_{config['gauges']}", 
        localrule: False
        group: f"evaluate_L{_level}"   #TODO: add this group to the config
        output: 
            performance = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "performance.nc"),
            eval_done = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "evaluate.done")
        threads: 1
        resources: 
            time = "00:30:00",
            mem_mb = 8000
        script: 
            """src/calib/evaluate_per_run.py""" # TODO: modify this script to submit multiple grouped jobs


    rule:
        name: f"combine_performance_L{_level}"
        input:
            done = expand(Path(calib_dir, f"level{_level}", '{params}', "evaluate.done"), params=paramspace.instance_patterns),
            performance_files=expand(Path(calib_dir, f"level{_level}", "{params}", "performance.nc"), params=paramspace.instance_patterns)
        output: 
            performance = Path(calib_dir, f"level{_level}", "performance.nc"),
            best_params = Path(calib_dir, f"level{_level}", "best_params.csv") #defaults to best 10
        localrule: True
        threads: 4
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
        localrule: True
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
        performance = glob.glob(str(Path(calib_dir, "level*", "performance.nc"))) #expand(Path(calib_dir, "{level}", "performance.nc"), level=list(graph.keys()))
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