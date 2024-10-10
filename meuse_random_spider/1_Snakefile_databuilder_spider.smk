# Subworkflow databuilder
# This subworkflow is used to build the data for the calibration of the wflow model

#@michaelohanrahan
#@jingdeng

import json
import shutil
import os
import platform
from pathlib import Path
from snakemake.utils import Paramspace
import geopandas as gpd
from src.calib.create_set_params import create_set
from glob import glob

#=========================================================
# Spider specific stuff
DRIVE = "/project"
PLATFORM = "Linux"

## Some preparatory actions
# Ensure the model directory is there
configfile: 'config/calib_spider.yml'

basin = config["basin"]
base_dir = Path(str(f'{DRIVE}/{config["base_dir"]}')).as_posix()                         # p: or /p/ ... / RWSOS_Calibration / basin
os.chdir(os.path.join(base_dir, basin))  # Set the working directory to the base directory
#setting as workdir with snakemake will look for scripts in the workdir not the executing directory
# workdir: base_dir

source_dir = Path(base_dir, basin, config["source_dir"])    # data/1-external
inter_dir = Path(base_dir, basin, config["inter_dir"])      # data/2-intermediate
input_dir = Path(base_dir, basin, config["input_dir"])      # data/3-input
out_dir = Path(base_dir, basin, config["output_dir"])       # data/4-output
vis_dir = Path(base_dir, basin, config["vis_dir"])          # data/5-visualisation

for dir in [source_dir, inter_dir, out_dir, vis_dir, input_dir]:
    os.makedirs(dir, exist_ok=True)

#Expected_levels
expected_levels = 6

# for dir in [source_dir, inter_dir, out_dir, vis_dir, input_dir]:
#     if not dir.exists():
#         os.makedirs(dir, exist_ok=True)

gauges = config["gauges"]

#=========================================================
#== Actions more focused on the calibration itself
# Load the dependency graph and sort into level
# Load parameter recipe as a parameter space
# ['KsatHorFrac', 'RootingDepth', 'SoilThickness'], ['set', 'add', 'mult'], df of params
rule all:
    input:  Path(inter_dir, "addksathorfrac", "staticmaps", "staticmaps.nc"),
            Path(inter_dir, f'{config["gauges"]}_levels_graph.json'),
            Path(inter_dir, "gaugesadded", f"wflow_sbm_addgauges.toml"),
            Path(inter_dir, "addflp", "wflow_sbm_addflp.toml"),
            Path(inter_dir, "addflp", "staticmaps", "staticmaps.nc"),
            Path(input_dir, "wflow_sbm.toml"),
            directory(expand(Path(inter_dir, "calib_data", "level{level}"), level=range(0, expected_levels))),
            Path(input_dir, "staticmaps", "staticmaps.nc"),
            Path(input_dir, "staticgeoms", f"subcatch_{config['gauges']}.geojson")
    default_target: True

'''
create the gaugemap from your preferred geojson file of gauges... best to make this with some manual checks to ensure the right gauges are being used
1: src/pre/assess_discharge_data.py to perform health checks and gather all gauges to be used (assuming all data is collected in a datacatalog)
2: scripts\pixi_run_gaugemap.bat will build a gaugemap from the output geojson and the dependency graph
 -- you can recursively modify the ignorelist in the interim folder to remove gauges to simplify the graph
    #TODO: automatically optimise the graph
3: src/pre/create_discharge_data.py to to create the discharge dataset from the gauges geojson
adding default target to false to encourage snakemake to run these rules before defining the wildcards
'''


rule create_gaugemap:
    input:
        gauges = Path(inter_dir, "QGIS", "to_wflow", f"{config['gauges']}_gauges_all.geojson")
    output:
        gridfile = Path(inter_dir, "gaugesadded", "staticmaps", "staticmaps.nc"),
        subcatch = Path(inter_dir, f"gaugesadded/staticgeoms/subcatch_{config['gauges']}.geojson"),
        gauges = Path(inter_dir, f"gaugesadded/staticgeoms/gauges_{config['gauges']}.geojson"),
        toml = Path(inter_dir, "gaugesadded", "wflow_sbm_addgauges.toml")
    params:
        cwd=Path(os.getcwd()).as_posix(),
        config_root=Path(source_dir),
        new_root=Path(inter_dir, "gaugesadded"),
        mode="w+",
        ignorelist=Path(inter_dir, "ignore_list.txt"),
        basename=config["gauges"],
        index_col="wflow_id",
        max_dist=1000,
        crs="EPSG:4326",
        config_old="wflow_sbm_template.toml",
        config_new="wflow_sbm_addgauges.toml"
    shell:
        """
        pixi run python src/pre/update_gauges.py {params.cwd} {params.config_root} {input.gauges} \
            --new_root "{params.new_root}" \
            --mode "{params.mode}" \
            --basename "{params.basename}" \
            --index_col "{params.index_col}" \
            --ignore_list "{params.ignorelist}" \
            --snap_to_river True \
            --max_dist {params.max_dist} \
            --derive_subcatch True \
            --crs "{params.crs}" \
            --config_old "{params.config_old}" \
            --config_new "{params.config_new}"
        """

"""
create the dependency graph from the gaugesadded gridfile, in this case is based off of the
hourly observations in the basin. This graph defines the elements for calculation and their
subsequent dependencies at each level, saving them in a json file for later use.
#TODO: establish the elements list for the expand function later
"""

rule create_graph:
    input:
        gridfile = rules.create_gaugemap.output.gridfile,
        toml = rules.create_gaugemap.output.toml,
        gauges = rules.create_gaugemap.output.gauges,
    params:
        gaugeset = config["gauges"]
    output:
        graph = Path(inter_dir, f'{config["gauges"]}_levels_graph.json')
    script:
        """src/graph/create_dependency_graph.py"""

rule create_folders:
    input:
        graph = rules.create_graph.output.graph
    # params:
    output:
        directory(expand(Path(inter_dir, "calib_data", "level"+"{level}"), level=range(0, expected_levels)))
    run:
        calib_dir = Path(inter_dir, "calib_data")
        os.makedirs(calib_dir, exist_ok=True)
        with open(input.graph, "r") as f:
            graph = json.load(f)

        # Ensure all the level directories are there
        # for _l in graph.keys():
        levels= []
        for level in graph.keys():
            levels.append(level)
            _p = Path(calib_dir, level)
            if not _p.exists():
                os.makedirs(_p, exist_ok=True)

#requires p drive connection for deltares data catalog
rule ksat_setup:
    input:
        gridfile_in=Path(inter_dir, "gaugesadded", "staticmaps", "staticmaps.nc"),
    output:
        config_fn_out=Path(inter_dir, "addksathorfrac", "wflow_sbm_addksathorfrac.toml"),
        gridfile_out = Path(inter_dir, "addksathorfrac", "staticmaps", "staticmaps.nc"),
        out_geoms = directory(Path(inter_dir, "addksathorfrac", "staticgeoms")),
        gauges = touch(Path(inter_dir, "addksathorfrac", "staticgeoms", f"gauges_{config['gauges']}.geojson")),
    params:
        drive=DRIVE,
        mod_new_root=Path(inter_dir, "addksathorfrac"),
        mod_root=Path(inter_dir, "gaugesadded"),
        config_fn_in="wflow_sbm_addgauges.toml",
        var=config['ksathorfrac_map'],
        config_fn_out="wflow_sbm_addksathorfrac.toml"
    shell:
        """
        pixi run python src/pre/setup_ksathorfrac.py \
            --DRIVE {params.drive} \
            --mod_root {params.mod_root} \
            --mod_new_root {params.mod_new_root} \
            --config_fn_in {params.config_fn_in} \
            --var {params.var} \
            --config_fn_out {params.config_fn_out}
        """

rule flp_setup:
    input:
        gridfile_in=Path(inter_dir, "addksathorfrac", "staticmaps", "staticmaps.nc"),
        config_fn_in=Path(inter_dir, "addksathorfrac", "wflow_sbm_addksathorfrac.toml"),
        geoms = Path(inter_dir, "addksathorfrac", "staticgeoms", f"gauges_{config['gauges']}.geojson")
    output:
        config_fn_out=Path(inter_dir, "addflp", "wflow_sbm_addflp.toml"),
        gridfile_out = Path(inter_dir, "addflp", "staticmaps", "staticmaps.nc"),
        geoms_out = directory(Path(inter_dir, "addflp", "staticgeoms"))
    params:
        config_fn_in="wflow_sbm_addksathorfrac.toml",
    script:
        """src/pre/insert_N_flp_var.py"""

rule new_model_to_input:
    input:
        config_fn=rules.flp_setup.output.config_fn_out,
        staticmaps=rules.flp_setup.output.gridfile_out
    params:
        new_name_toml = Path(input_dir, "wflow_sbm.toml"),
        new_name_nc = Path(input_dir, "staticmaps", "staticmaps.nc")
    output:
        config_run=Path(input_dir, "wflow_sbm.toml"),
        staticmaps_run=Path(input_dir,"staticmaps", "staticmaps.nc")
    run:
        shutil.copy(input.config_fn, params.new_name_toml)
        shutil.copy(input.staticmaps, params.new_name_nc)

rule staticgeoms_to_input:
    input:
        staticmaps=rules.create_gaugemap.output.gridfile
    params:
        sg_folder = Path(source_dir, "staticgeoms"),
        new_gauges = Path(inter_dir, "gaugesadded", "staticgeoms", f"gauges_{config['gauges']}.geojson"),
        subcatch = Path(inter_dir, "gaugesadded", "staticgeoms", f"subcatch_{config['gauges']}.geojson")
    output:
        sg_sc=Path(input_dir, "staticgeoms", f"subcatch_{config['gauges']}.geojson"),
        sg_g=Path(input_dir, "staticgeoms", f"basins.geojson")
    run:
        # Ensure the destination directory exists
        if os.path.exists(Path(input_dir, "staticgeoms")):
            shutil.rmtree(Path(input_dir, "staticgeoms"))
        os.makedirs(Path(input_dir, "staticgeoms"), exist_ok=True)

        # Copy all files from the source directory to the destination directory
        for item in os.listdir(params.sg_folder):
            s = Path(params.sg_folder, item)
            d = Path(input_dir, "staticgeoms", item)
            if s.is_file():
                shutil.copy2(s, d)

        # Then copy the new gauges
        shutil.copy(params.new_gauges, output.sg_g)
        shutil.copy(params.subcatch, output.sg_sc)

rule lakes_to_input:
    input:
        lakes = expand(Path(source_dir, "staticmaps", "lake_hq_"+"{hq}"+".geojson"), hq=[1, 2, 3])
    output:
        lakes = expand(Path(input_dir, "staticmaps", "lake_hq_"+"{hq}"+".geojson"), hq=[1, 2, 3])
    run:
        os.makedirs(Path(input_dir, "staticmaps"), exist_ok=True)
        shutil.copy(input.lakes, output.lakes)