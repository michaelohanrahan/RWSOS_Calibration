import json
import shutil
from pathlib import Path

from snakemake.utils import Paramspace

from src.calib.create_set_params import create_set
from src.calib.dependency_graph import sort_graph


## Some preparatory actions
# Ensure the model directory is there
base_dir = Path(config["base_dir"])
root_dir = Path(base_dir, config["output_dir"])
source_dir = Path(base_dir, config["source_dir"])
if not root_dir.exists():
    root_dir.mkdir()

# Copy the needed files to be filled in
cfg_template = Path(root_dir, config["base_config"])
staticmaps = Path(root_dir, "staticmaps.nc")
shutil.copy2(
    Path(source_dir, config["base_config"]),
    cfg_template,
)
if not staticmaps.exists():
    shutil.copy2(
        Path(source_dir, config["base_staticmaps"]),
        staticmaps,
    )

# Ensure the folder is there for all the calibration data and settings
calib_dir = Path(root_dir, "calib_data")
if not calib_dir.exists():
    calib_dir.mkdir()


## Actions more focussed on the calibration itself
# Load the dependency graph and sort into levels
graph = sort_graph(
    Path(config["dependency_graph"])
)
# Ensure all the level directories are there
for _l in graph.keys():
    _p = Path(calib_dir, _l)
    if not _p.exists():
        _p.mkdir()

# Get the last index of the keys
last_level = list(graph.keys())[-1]

# Get all the elements for the visualization
elements = []
for _l in graph.values():
    elements += _l["elements"]

# Load parameter recipe as a parameter space
lnames, methods, ds = create_set(config["calib_recipe"])
paramspace = Paramspace(ds, filename_params="*", param_sep="", filename_sep="_")

############################
# DOING the snakey!
############################
# Define the 
rule all:
    input: expand(Path(root_dir, config["result_dir"], "figures", "hydro_{gauge}.png"), gauge=elements)

rule config:
    input: lambda wildcards: expand(Path(calib_dir, "{dep}", "done.txt").as_posix(), dep=graph[wildcards.item]["deps"]) if graph[wildcards.item]["deps"] else []
    params: 
        cfg_template = cfg_template,
        starttime = config["starttime"],
        endtime = config["endtime"],
        timestep = config["timestep"],
        forcing_path = Path(source_dir, config["source_forcing_data"])
    output: expand(Path(calib_dir, "{item}", "{params}", config["wflow_cfg_name"]), item=f"{{item}}", params=paramspace.wildcard_pattern)
    script:
        """src/calib/set_config.py"""

rule create_params:
    input: Path(calib_dir, "{item}", paramspace.wildcard_pattern, config["wflow_cfg_name"])
    params:
        dataset = staticmaps,
        params = paramspace.instance,
        params_lname = lnames,
        params_method = methods,
        level = "{item}",
        graph = graph,
        sub_catch = Path(source_dir, config["subcatch"])
    output: 
        staticmaps = Path(calib_dir, "{item}", paramspace.wildcard_pattern, "staticmaps.nc")
    script: 
        """src/calib/set_calib_params.py"""

rule wflow:
    input: 
        cfg = Path(calib_dir, "{item}", paramspace.wildcard_pattern, config["wflow_cfg_name"]),
        staticmaps = Path(calib_dir, "{item}", paramspace.wildcard_pattern, "staticmaps.nc")
    output: Path(calib_dir, "{item}", paramspace.wildcard_pattern, "run_default", "output_scalar.nc")
    shell: 
        f"""julia --project="{config['wflow_project_dir']}" -t {config['wflow_threads']} -e "using Wflow; Wflow.run()" {{input.cfg}}"""

rule evaluate:
    input: expand(Path(calib_dir, "{item}", "{params}", "run_default", "output_scalar.nc"), item=f"{{item}}", params=paramspace.instance_patterns)
    params: 
        observed_data = config["observed_data"],
        level = "{item}",
        graph = graph,
        params = ds.to_dict(orient="records"), 
        starttime = config["eval_starttime"],
        endtime = config["eval_endtime"],
        metrics = config["metrics"],
        weights = config["weights"]
    output: 
        best_params = Path(calib_dir, "{item}", "best_params.csv"),
        performance = Path(calib_dir, "{item}", "performance.nc")
    script: 
        """src/calib/evaluate_params.py"""

rule set_params:
    input: 
        best_params = Path(calib_dir, "{item}", "best_params.csv")
    params:
        staticmaps = staticmaps,
        sub_catch = Path(source_dir, config["subcatch"]),
        params_lname = lnames,
        params_method = methods
    output: 
        done = Path(calib_dir, "{item}", "done.txt")
    script:
        """src/calib/set_eval_params.py"""

rule prep_final_stage:
    input: 
        done = Path(calib_dir, last_level, "done.txt"),
        performance = expand(Path(calib_dir, "{level}", "performance.nc"), level=list(graph.keys()))
    params:
        cfg_template = cfg_template,
        cfg_args = [config["starttime"], config["endtime"], config["timestep"], Path(source_dir, config["source_forcing_data"])],
        staticmaps = staticmaps
    output: 
        cfg = Path(root_dir, config["result_dir"], config["wflow_cfg_name"]),
        performance = Path(root_dir, config["result_dir"], "performance.nc"),
        staticmaps = Path(root_dir, config["result_dir"], "staticmaps.nc")
    script:
        """src/calib/prep_final_stage.py"""

rule run_final_model:
    input:
        cfg = Path(root_dir, config["result_dir"], config["wflow_cfg_name"]),
        staticmaps = Path(root_dir, config["result_dir"], "staticmaps.nc")
    output: Path(root_dir, config["result_dir"], "run_default", "output_scalar.nc")
    shell:
        f"""julia --project="{config['wflow_project_dir']}" -t {config['model_threads']} -e "using Wflow; Wflow.run()" {{input.cfg}}"""

rule visualize:
    input: 
        scalar = Path(root_dir, config["result_dir"], "run_default", "output_scalar.nc"),
        performance = Path(root_dir, config["result_dir"], "performance.nc")
    params:
        observed_data = config["observed_data"],
        gauges = elements,
        starttime = config["eval_starttime"],
        endtime = config["eval_endtime"],
        period_startdate = config["hydro_period_startdate"],
        period_length = config["hydro_period_length"],
        period_unit = config["hydro_period_unit"],
        output_dir = Path(root_dir, config["result_dir"], "figures")
    output:
        figures = expand(Path(root_dir, config["result_dir"], "figures", "hydro_{gauge}.png"), gauge=elements)
    script:
        """src/post/plot_final_model.py"""