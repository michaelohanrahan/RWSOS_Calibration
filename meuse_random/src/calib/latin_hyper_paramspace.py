from scipy.stats import qmc
import json
<<<<<<< HEAD
from .setuplog import setup_logging
=======
from setuplog import setup_logging
>>>>>>> 50f0b2a433687946660e04c31a5fc4d5b3a8b94a
import traceback
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def generate_samples(
                    #  l, #logger
                     LEVEL:int,
                     RECIPE:str|Path, 
                     N_SAMPLES:int,
                    #  OUT:str|Path,
                     OPTIM=None,):
    
    #calibration recipe
    with open(RECIPE, 'rb') as f:
        data = json.load(f)

    #Seed per level to acheive reproducibility
    seed = np.random.default_rng(LEVEL)
    
    l_names = []
    s_names = []
    methods = []
    l_bounds_list = []
    u_bounds_list = []

    for k, v in data.items():
        l_names.append(k)
        s_name = v['short_name']
        methods.append(v['method'])
        s_names.append(s_name)
        l_bounds_list.append(min(v['values']))
        u_bounds_list.append(max(v['values']))


    #initialize the Latin Hypercube sampler
    lhd = qmc.LatinHypercube(d=len(s_names),
                             optimization=OPTIM, 
                             seed=seed
                        )
    #Create a sampler
    samples = lhd.random(n=N_SAMPLES)

    #Scale the samples to the bounds
    sample_scaled = qmc.scale(samples, l_bounds_list, u_bounds_list)
    
    #Create a dataframe
    params_df= pd.DataFrame(sample_scaled, columns=s_names, index=[LEVEL]*N_SAMPLES)
    
    #Split any columns with comma separation (co-scaling)
    cols = []
    for s_name in s_names:
        if ',' in s_name:
            col1, col2 = s_name.split(',')
            cols.append(col1)
            cols.append(col2)
        else:
            cols.append(s_name)
            
    # Handle columns with comma separation
    for s_name in s_names:
        # print(s_name)
        if ',' in s_name:
            col1, col2 = s_name.split(',')
            params_df[col1] = params_df[s_name]
            params_df[col2] = params_df[s_name]
            params_df.drop(columns=s_name, inplace=True)
    
    params_df = params_df[cols]  # Reorder columns to match the original order
    params_df.index.name = 'level'
    
    return l_names, methods, params_df.round(2)


def create_set_all_levels(last_level,
                          RECIPE,
                          N_SAMPLES,
                          OPTIM,
                          ):
    
    _df_list = []
    for LEVEL in range(0, last_level+1):
        l_names, methods, params_df = generate_samples(
                                    LEVEL=LEVEL, 
                                    RECIPE=RECIPE, 
                                    N_SAMPLES=N_SAMPLES, 
                                    OPTIM=OPTIM,                            
        )
        _df_list.append(params_df)
        
    all_params_df = pd.concat(_df_list, axis=0)
    
    return l_names, methods, all_params_df



if __name__ == "__main__":
    
    # test dev
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\UNREAL_TEST_DATA')
    
    # find last level from the final level directory
    import glob
    levels = glob.glob(str(Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim','calib_data', "level*")))
    levels_ints = [int(level.split("level")[-1]) for level in levels]
    last_level = int(levels[-1].split("level")[-1])
    
    # set up input vars
    l=setup_logging(work_dir, 'latin_hyper_paramspace.log')
    RECIPE = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\meuse\config\calib_recipe.json')
    N_SAMPLES = 1000
    OPTIM = 'random-cd'
    # OUT = work_dir / 'LHS_df.csv'
    
    l_names, methods, all_params_df = create_set_all_levels(
                          last_level,
                          RECIPE,
                          N_SAMPLES,
                          OPTIM,
                          )
    
    
    
    l=setup_logging('data/0-log', 'latin_hyper_paramspace.log')
    try:
        LEVEL = 0
        RECIPE = '../../config/calib_recipe.json'
        N_SAMPLES = 1000
        OPTIM = 'random-cd'
        OUT=Path('../../data/2-interim/calib_data/level0/params.csv')
        l_names, methods, params_df = generate_samples(l,
                                    LEVEL=LEVEL, 
                                    RECIPE=RECIPE, 
                                    N_SAMPLES=N_SAMPLES, 
                                    OPTIM=OPTIM,
                                    OUT=OUT)
        
        plt.figure(figsize=(10, 6))
        params_df.boxplot()
        plt.title('Boxplot of the parameters')
        plt.savefig('../../data/2-interim/calib_data/level0/params_boxplot.png')
        
        #3d scatter for the first three parameters
        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(params_df.iloc[:,0], params_df.iloc[:,1], params_df.iloc[:,2])
        ax.set_xlabel(params_df.columns[0])
        ax.set_ylabel(params_df.columns[1])
        ax.set_zlabel(params_df.columns[2])
        plt.title("Scatter of the first three parameters")
        plt.savefig('../../data/2-interim/calib_data/level0/params_3dscatter.png')
        
        #histogram of the first parameter
        fig = plt.figure(figsize=(10, 6))
        plt.hist(params_df.iloc[:,0], bins=20, color='blue', edgecolor='black')
        plt.title('Histogram of the first parameter')
        plt.savefig('../../data/2-interim/calib_data/level0/params_hist.png')
        
    except Exception as e:
        l.error(f'An error occurred: {e}')
        l.error(traceback.format_exc())
        raise e
        
    