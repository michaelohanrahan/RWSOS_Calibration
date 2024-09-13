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
import ast

def fill_pred_params_recursive(pred, 
                     chosen_param, 
                     pred_level, 
                     graph_pred, 
                     random_params_all, 
                     random_params):
    """
    pred: float, the gauge id of the predecessor currently being filled
    chosen_param: dict, the paramset of the current gauge
    pred_level: str, the level of the predecessor
    graph_pred: dict, the graph of the predecessors
    random_params_all: dict, all the random params for each level
    random_params: pd.DataFrame, the random params for the current level
    """
    #TODO: handle a number of predecessors
    #THEN: RETURN THE pre_pre
    # the level is not as minimized, for this paramset we need the old choices added
    random_params_level = random_params_all[pred_level]
    _mask = pd.Series([True] * len(random_params_level), index=random_params_level.index)
    
    for key, value in chosen_param.items():
        _mask = _mask & (random_params_level[key] == value)
    
    #row matching the params of the current gauge, containing the random params for the pred
    params_pred = random_params_level[_mask]
    
    #now we need to fill in the random params for the pred's pred and so on
    pred_upgauge_ids = graph_pred[str(pred)]['_pre']  # list of float
    
    remaining_upgauge = []
    done_gauge = []
    for pred_up in pred_upgauge_ids:
        
        #The ones we are working with are filled here
        random_params[pred_up] = params_pred.loc[pred_up]
        pred_level = graph_pred[str(pred_up)]['level']
        
        fill_pred_params_recursive(pred_up, 
                                    chosen_param, 
                                    pred_level, 
                                    graph_pred, 
                                    random_params_all, 
                                    random_params)
        
    return random_params, done_gauge, remaining_upgauge


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
    """
    Generating random parameter combinations for the immediate upstream gauges.
    Random integer values select from the top n sets from the predecessors.
    Using predecessors because the current gauge could have a variety of levels as immediate upstream
    for example: main basin outflow (max level) could have a predecessor at any preceding level. 
    
    Best params previous is a concatenation of all previous levels, index by level and gauge.
       -- This way we can access the best of any level 
       -- The upgauge random params are generated as new columns for in the paramspace df 
       -- resulting in a random_params df that will successfully retain all the correct values
    
    This is sound for L0 to L1, but L1, L2 ?? The information persists, but there is no master combination sheet. 
    
    I think this qualifies as a TODO:
    - we should envision the random params as a df with a column for each predecessors predecessor. 
    - and then we fill in the values for each level. 
    """
    if level != 'level0':
        # load best_params_previous
        best_params_previous = pd.read_csv(best_params_previous, index_col=['level', 'gauge'])
        level_int = int(level.split('level')[-1])
        
        # e.g. @ level3: "level1": df(level1/random_params.csv)}, {"level2": df(level2/random_params.csv)
        if level_int > 1:
            bp_format = Path(best_params_previous).parent.parent / "level{}" / "random_params.csv"
            random_params_all = {
                f"level{_l}":pd.read_csv(bp_format.format(_l)) for _l in range(1,level_int)
            }
        
        
        random_params = pd.read_csv(params_df)
        
        # current level gauges
        gauges = graph[level]['elements']  # list of float
        for g in gauges:
            # select direct upstream gauges (graph_pred)
            pred_upgauge_ids = graph_pred[str(g)]['_pre']  # list of float
            
            # search in the best_params_previous for this direct upstream
            # set predecessors randomly
            for pred in pred_upgauge_ids:  # pred: float
                pred_level = graph_pred[str(pred)]['level']
                
                # how many best and get random one
                n_Top = [int(i.split('_')[-1]) for i in best_params_previous.columns if 'Top' in i]
                Top_max = max(n_Top)
                random_col = random.randint(1, Top_max)
                
                # select the random one and fill in random_params dataframe
                chosen_param = best_params_previous.loc[(pred_level, int(pred)), f"Top_{random_col}"]
                
                #this will give the paramset for this gauge
                random_params[pred] = chosen_param
                
                #now we need to fill in the random params for the pred's pred and so on
                dep_upgauge_ids = graph_pred[str(pred)]['_pre']  # list of float
                
                #still remaining 
                remaining_upgauge = [graph_node[str(i)]['_deps'] for i in dep_upgauge_ids]
                remaining_upgauge = [item for sublist in remaining_upgauge for item in sublist]
                
                for d_pred in dep_upgauge_ids:
                    #Setting the determined params
                    d_pred_level = graph_pred[str(d_pred)]['level']
                    random_params, done_gauge, remaining_upgauge = fill_pred_params_recursive(
                        d_pred, 
                        d_pred_level, 
                        graph_pred, 
                        random_params_all, 
                        random_params
                        )
                    
                      
        
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
                
                lnames, methods, all_level_df = create_set_all_levels(last_level=5, RECIPE="config/LHS_calib_recipe.json", N_SAMPLES=10, OPTIM='random-cd')
                graph = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_levels_graph.json"))
                graph_pred = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_pred_graph.json"))
                graph_node = json.load(open("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/Hall_nodes_graph.json"))
                level = 'level1'
                best_params_previous = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level0/best_params.csv"
                params_df = "/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data/level1/paramspace.csv"
                out = Path('data/2-interim','calib_data', level, 'random_params.csv')
                
                
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


