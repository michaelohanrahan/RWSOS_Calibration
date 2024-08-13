import json
import pickle as pk
import pandas as pd
from pathlib import Path
import random
import xarray as xr
from setuplog import setup_logging
import traceback
import os 

def main(
    l,
    level: str,
    best_params: pd.DataFrame,
    params_df: pd.DataFrame,
    graph: dict,
):
    #TODO: as change to LHS, params_df should be the param space for the next level
    random_df = params_df
    gauges = graph[level]['elements']
    level_int = int(level.split('level')[-1])
    random_df['level'] = level_int
    
    #how many best
    n_Top = [int(i.split('_')[-1]) for i in best_params.columns if 'Top' in i]
    Top_max = max(n_Top)
    
    for gauge in gauges:
        random_df[gauge] = None
        random_df[gauge] = random_df[gauge].apply(lambda x: best_params.loc[gauge, f"Top_{random.randint(1, Top_max)}"])
    
    out_dir = Path(f"data", "2-interim", "calib_data", level)
    os.makedirs(out_dir, exist_ok=True)
    random_out = Path(out_dir, "random_params.csv")
    
    random_df.to_csv(random_out, index=False)
    l.info(f"Random parameters saved to {random_out}")
    with open(Path(out_dir, "done.txt"), "w") as f:
        f.write("done")
    l.info(f"Done flag saved to {out_dir/'done.txt'}")

if __name__ == "__main__":
    # Set up logger
    l = setup_logging('data/0-log', f'random_params.log')
    try:
        if 'snakemake' in globals():
            snakemake = globals()["snakemake"]
            best_params = snakemake.input.best_params #"best_n_params"
            level = snakemake.params.level #"leveln"
            params_df = snakemake.params.params #"params_df"
            graph = snakemake.params.graph
        else:
            try:
                cwd = Path(r"c:\Users\deng_jg\work\05wflowRWS")
                os.chdir(cwd)
                data_dir = Path('UNREAL_TEST_DATA')
                p = data_dir / 'staticmaps/staticmaps.nc'
                
                with open(data_dir/'create_set_params.pkl', 'rb') as f:
                    d = pk.load(f)
                
                params_lname = d['lnames']
                params_method = d['methods']
                params_df = d['ds']
                graph = json.load(open(Path(data_dir, 'Hall_levels_graph.json')))
                level = 'level0'
                gauges = graph[level]['elements']
                data = {
                    'gauges': list(gauges),
                    **{'Top_{}'.format(i): str(dict(params_df.iloc[random.randint(0, len(params_df)-1)])) for i in range(10+1)},
                }
                best_params = pd.DataFrame(data)
                best_params = best_params.set_index('gauges')
                
                
            except Exception as e:
                l.error("An error occurred while setting default arguments", exc_info=True)
                raise e

        l.info(f'cwd: {os.getcwd()}')
        l.info(f'logger: {l}')
        main(
            l,
            level=level,
            best_params=best_params,
            params_df=params_df,
            graph=graph,
        )
        
    except Exception as e:
        l.error(f"An error occurred: {e}")
        l.error(traceback.format_exc())
        raise e


