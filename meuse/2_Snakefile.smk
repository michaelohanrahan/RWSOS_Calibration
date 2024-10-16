
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

#Base directory
basin = config["basin"]        # meuse
base_dir = config["base_dir"]  # RWSOS_Calibration
base_dir = f"{DRIVE}/{Path(base_dir).as_posix()}" # p: or /p/ ... / RWSOS_Calibration / basin
workdir: str(Path(base_dir, basin).as_posix())

# Workflow for the calibration of the wflow model
# This workflow is used to calibrate the wflow model using the data built in the databuilder subworkflow
#@michaelohanrahan
#@jingdeng

import json
import shutil
import os
from snakemake.utils import Paramspace
import geopandas as gpd
from src.calib.create_set_params import create_set
import glob

#=========================================================
log_dir = Path(base_dir, basin, config["log_dir"]).as_posix()          # data/0-log
source_dir = Path(base_dir, basin, config["source_dir"]).as_posix()    # data/1-external
inter_dir = Path(base_dir, basin, config["inter_dir"]).as_posix()      # data/2-intermediate
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

#paramspace
lnames, methods, df = create_set(config["calib_recipe"])
paramspace = Paramspace(df)

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

#graph
graph = json.load(open(Path(inter_dir, f'{config["gauges"]}_levels_graph.json')))

#soilthicknesses
with open(config["calib_recipe"]) as recipe:
    recipe = json.load(recipe)
    # print(recipe)
    ST_values = list(recipe["SoilThickness_manual_cal"]["values"])
    ST_str = [str(ST).replace('.', '') for ST in ST_values]
    ST_dict = {ST_str[i]: ST_values[i] for i in range(len(ST_values))}

############################
# DOING the snakey!
############################
# Define the main rule all that expects all the visualization files to have been created for successfull completion
# rule all:
# Define the main rule all that expects all the visualization files to have been created for successful completion
#TODO: Figure out dynamic resources from the slurm availability
'''
:: rule all ::
Defines the overall default target
We expect upon successfull completion that all gauges will have been visualized
'''
rule all:
    input: 
        expand(Path(vis_dir, "hydro_gauge", "hydro_{gauge}.png"), gauge=elements),
        expand(Path(input_dir, "instates", "instate_level{level}_ST{_t_str}.nc"), level=range(0, last_level+1), _t_str=ST_str),
        expand(Path(calib_dir, "level{nlevel}", "done.txt"), nlevel=range(-1, last_level+1)),
    default_target: True

#THIS RULE ALLOWS RECURSIVE DEPENDENCY
rule init_done:
    params:
        d = Path(calib_dir, "level-1")
    output: 
        p = Path(calib_dir, "level-1", "done.txt")
    localrule: True
    run:
        os.makedirs(params.d, exist_ok=True)
        with open(output.p, "w") as f:
            f.write("done")

#Had to switch to a looping approach to get instates per level
#Otherwise they would not wait for the previous level to finish
for _level in range(5, last_level+1):
    '''
    :: initial instates ::
    This rule isnt waiting for the output of a previous rule, so it is run first and only once.
    That is why we hardcode the level. The subsequent rules will have to wait for an input to run for instates.
    #Try running this for 0 and then later dynamically make TOMLs for each level that react to the output of the main runs
    '''
    rule: 
        name: f"initial_instate_tomls_L{_level}"
        input:
            config_fn = cfg_template, #the input template
            prev = Path(calib_dir, f"level{_level-1}", "done.txt") #the output of the previous level
        params:
            ST_values = ST_values,
            root = Path(input_dir).as_posix(),
            ST_key= config["soilthickness_map"],
            level = _level,
            starttime = config["instate_starttime"],
            endtime = config["instate_endtime"],
        output: #TOML expects: f"staticmaps_L{level}_ST{soilthickness}.nc"
            ST_instates_tomls = expand(Path(input_dir, "instates", "wflow_sbm_getinstate_level"+f"{_level}"+"_ST{_t_str}.toml"), _t_str=ST_str)
        localrule: True
        script:
            """src/calib/create_instate_toml.py"""
    '''
    Now in the looping structure we can wait for the best parameters upstream before we run the instates
    the downside curently is the staticmaps in the input is consistently overwritten per level... 
    '''
    rule:
        name: f"ST_instate_staticmaps_L{_level}"
        input: 
            in_tomls = Path(input_dir, "instates", "wflow_sbm_getinstate_level"+f"{_level}"+"_ST{_t_str}.toml")
        params:
            base_grid = staticmaps,
            level = _level,
            ST_str = "{_t_str}",
            ST_key = config["soilthickness_map"],
            graph = graph,
            sub_catch = subcatch,
            lakes_in = lakes
        localrule: True
        output:
            ST_grids = Path(input_dir, "instates", "staticmaps_L"+f"{_level}"+"_ST{_t_str}.nc")
        script:
            """src/calib/set_instate_ST.py"""
    
    '''
    :: Run instates::
    We wait for the output of the initial state defining tomls and then we run the instates
    This will hopefully pick up on level 0, and then later on the other levels, hardcoding 0 here eliminates recursion
    '''
    rule:
        name: f"run_instate_L{_level}"
        input:
            grid = Path(input_dir, "instates", "staticmaps_L"+f"{_level}"+"_ST{_t_str}.nc"),  
            cfg = Path(input_dir, "instates", "wflow_sbm_getinstate_level"+f"{_level}"+"_ST{_t_str}.toml"),
        params:
            project = Path(base_dir, "bin").as_posix(),
            threads = config["wflow_threads"]
        output:
            done = Path(input_dir, "instates", f"instate_level{_level}_ST"+"{_t_str}.nc")
        group: f"instate_L{_level}"
        localrule: False
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
    
    """
    config: This rule modifies a blueprint configuration file for a specific time period and forcing path. 
    It takes the blueprint configuration file, start time, end time, time step, and forcing path as parameters. 
    The output is a set of modified configuration files.
    #DONE: in each config that accesses a certain level & soilthickness the instate reference should match
        If I can access the level and soilthickness from the {item} and {params} then we can point to the correct instate for each run
    """
    rule:
        name: f"config_L{_level}"
        # input: lambda wildcards: expand(Path(calib_dir, "{dep}", "done.txt").as_posix(), dep=graph[wildcards.item]["deps"]) if graph[wildcards.item]["deps"] else []
        input: expand(Path(input_dir, "instates", f"instate_level{_level}"+"_ST{_t_str}.nc"), _t_str=ST_str) #, level=f"{{item}}", thickness=[str(ST).replace('.', '') for ST in ST_values])
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
    create_params:  This rule creates parameter sets for the Wflow model. 
                    It takes a configuration file, a dataset of static maps, a parameter space instance, 
                    a list of parameter names, a list of parameter methods, a level, a graph, and a sub-catchment as parameters. 
                    The output is a set of static map files.
    #ADDED: lake_hqs to the output since i cannot see how to specify their specific location they are now alognside each grid file
    '''

    rule:
        name: f"create_params_L{_level}"
        input: Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, config["wflow_cfg_name"])
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
    wflow: This rule runs the hydrological model. It takes a configuration file and a grid (staticmap) file as input.
    #DONE: Wflow expects the appropriate instate
    #TODO: Only using the CAL time split...??? Shorter run time worse model...
    #TODO: modify the templace output to do output_scalar.nc
    '''

    rule:
        name: f"wflow_L{_level}"
        input: 
            instates=expand(Path(input_dir, "instates", f"instate_level{_level}"+"_ST"+"{_t_str}"+".nc"), _t_str=ST_str),
            cfg = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, config["wflow_cfg_name"]),
            staticmaps = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "staticmaps.nc"),
            lake_hqs = expand(Path(calib_dir, f"level{_level}", "{params}", "{lakes}"), params=paramspace.wildcard_pattern, lakes=lakefiles)
        params:
            project = Path(base_dir, "bin").as_posix(),
            # threads = config["wflow_threads"]
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
    Output:  an unstacked performance.nc file with metrics and weights for each parameter set, and a best_params.csv file with the best parameter set.
            out csv with the best parameter set for each gauge 
            
    '''

    rule:
        name: f"evaluate_L{_level}"
        input: 
            sim = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "output_scalar.nc"),
            done = Path(calib_dir, f"level{_level}", paramspace.wildcard_pattern, "run.done"),
        params: 
            observed_data = Path(source_dir, config["observed_data"]),
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
            performance = directory(Path(calib_dir, f"level{_level}", "performance.zarr")),
            best_params = Path(calib_dir, f"level{_level}", "best_params.csv") #defaults to best 10
        localrule: False
        threads: 4
        resources:
            time = "12:00:00",
            mem_mb = 32000
        script:
            "src/calib/combine_evaluated.py"

    '''
    This rule overwrites the staticmaps file with the best per level??
    #TODO: create a staticmaps per level? Not necessary if we are waiting for the done.txt
            -- not necessary but as a failsafe we create a copy of the staticmaps for each level
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
            done = Path(calib_dir, f"level{_level}", "done.txt")
        script: 
            """src/calib/set_eval_params.py"""

'''
Preparing the final stage: This rule prepares the final stage of the calibration process.
'''

rule prep_final_stage:
    input: 
        done = Path(calib_dir, "level"+f'{last_level}', "done.txt"),
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