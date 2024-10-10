import json
import pandas as pd
from pathlib import Path
import random
import xarray as xr
from setuplog import setup_logging
from latin_hyper_paramspace import create_set_all_levels
import traceback
import os 
from icecream import ic

def main(
    l,
    level: str,
    best_params_previous: Path | str,
    params_df: pd.DataFrame,  #DONE: as change to LHS, params_df should be the param space for this level
    graph: dict,  # Hall_levels_graph.json
    graph_pred: dict,  # Hall_pred_graph.json
    graph_node: dict,  # Hall_nodes_graph.json 
    out: Path | str,   #random_params.csv
):
    
    if level != 'level0':
        # load best_params_previous
        best_params_previous = pd.read_csv(best_params_previous, index_col=['level', 'gauge'], dtype={'gauge': int})
        random_params = pd.read_csv(params_df)
        
        # current level gauges
        gauges = graph[level]['elements']  # list of float
        
        for g in gauges:
            # select direct upstream gauges (graph_pred)
            pred_upgauge_ids = graph_pred[str(g)]['_pre']  # list of float
            
            # search in the best_params_previous for this direct upstream
            for pred in pred_upgauge_ids:  # pred: float
                pred_level = graph_pred[str(pred)]['level']  # the level of this direct upstream gauge
                upgauges_ids = graph_node[str(pred)]['_deps']  # list of further upstream gauges of this direct upstream gauge (float
                
                # Generate a list of 3000 random_col between 1 and Top_max
                n_Top = [int(i.split('_')[-1]) for i in best_params_previous.columns if 'Top' in i]
                Top_max = max(n_Top)
                
                random.seed(int(level.lstrip('level'))+int(pred))
                random_cols = [random.randint(1, Top_max) for _ in range(len(random_params))]
                
                # For each row in random_params, select a random param for pred and upgauge from best_params_previous
                for index in range(len(random_params)):
                    random_col = random_cols[index]
                    random_params.at[index, pred] = best_params_previous.loc[(pred_level, int(pred)), f"Top_{random_col}"] # locate the direct upstream gauge
                    for upgauge in upgauges_ids:  # upgauge: float. locate the further upstream gauge
                        random_params.at[index, upgauge] = best_params_previous.loc[(pred_level, int(upgauge)), f"Top_{random_col}"]
                        
        random_params.to_csv(out, index=False)
        l.info(f"Random parameters saved to {out}")
    else:
        l.info(f"Level 0 does not have upstream gauges. Random parameters not generated.")
        random_params = {"gauges": list(graph[level]['elements'])}
        df = pd.DataFrame(random_params)
        df.to_csv(out, index=False)
        l.info(f"Blank data saved to {out}")
        
        

if __name__ == "__main__":
    
    # Set up logger
    l = setup_logging('data/0-log', f'01-random_params.log')
    try:
        if 'snakemake' in globals():
            snakemake = globals()["snakemake"]
            best_params_previous = snakemake.input.best_params_previous #"best_n_params"
            level = snakemake.params.level #"leveln"
            params_df = snakemake.params.params_df #"params_df"
            graph = snakemake.params.graph
            graph_pred = snakemake.params.graph_pred
            graph_node = snakemake.params.graph_node
            out = snakemake.output.random_params
        else:
            try:
                
                # p = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/3-input/staticmaps/staticmaps.nc"
                
                # lnames, methods, all_level_df = create_set_all_levels(last_level=5, RECIPE="config/LHS_calib_recipe.json", N_SAMPLES=10, OPTIM='random-cd')
                graph = json.load(open(r"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_levels_graph.json"))
                graph_pred = json.load(open(r"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_pred_graph.json"))
                graph_node = json.load(open(r"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_nodes_graph.json"))
                level = 'level1'
                level_previous = 'level0'
                best_params_previous = rf"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random\data\2-interim\calib_spider\{level_previous}\best_params.csv"
                params_df = rf"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random\data\2-interim\calib_spider\{level}\paramspace.csv"
                out = Path(f"p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_spider/{level}/random_params_jing_seed2.csv")
                
                
            except Exception as e:
                l.error("An error occurred while setting default arguments", exc_info=True)
                raise e

        l.info(f'cwd: {os.getcwd()}')
        l.info(f'logger: {l}')
        
        main(
            l,
            level=level,
            best_params_previous=best_params_previous,
            params_df=params_df,
            graph=graph,
            graph_pred=graph_pred,
            graph_node=graph_node,
            out=out,
        )
        
        
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e


