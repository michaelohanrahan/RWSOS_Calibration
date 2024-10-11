from hydromt_wflow import WflowModel
from hydromt.data_catalog import DataCatalog
from hydromt.log import setuplog
import hydromt_wflow
import hydromt
import os 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import xarray as xr
from itertools import product


logger = setuplog("wflow", log_level=10,)
logger.info(f'Hydromt version: {hydromt_wflow.__version__}')
var1 = 'BRT_'
var2 = 'RF_'

mod_dict = {'base': r'p:\11209265-grade2023\wflow\RWSOS_Calibration\rhine\data\2-interim\ksat_test_base',
            'brt': r'p:\11209265-grade2023\wflow\RWSOS_Calibration\rhine\data\2-interim\ksat_test_brt',
            'rf': r'p:\11209265-grade2023\wflow\RWSOS_Calibration\rhine\data\2-interim\ksat_test_rf'}

ksat_dict = {'base': 'KsatHorFrac',
            'brt': 'ksathorfrac_BRT_250',
            'rf': 'ksathorfrac_RF_250'}

mod_dict = {}

for key, mod_root in mod_dict.items():
    
    mod = WflowModel(root=mod_root, 
                config_fn="wflow_sbm.toml",
                mode='r',
                logger=logger)
    
    mod.read_config()
    mod.read_grid()
    
    mod_dict[key] = mod

ksat_dict = {}	

for key, mod in mod_dict.items():
    ksat_dict[key] = mod.grid[ksat_dict[key]]


#%%
# Create a colormap that has transparency for values less than 0
cmap = plt.get_cmap('Reds')
cmap.set_under('none')

# Create a norm that maps all values less than 0 to the 'under' color
norm = colors.Normalize(vmin=0, vmax=5000)

maps = [ks1, ks2, ks3]
varnames = ['ksathorfrac_sub', 'ksathorfrac_BRT_250', 'ksathorfrac_RF_250']

# Plot each map
for map, var in zip(maps, varnames):
    plt.figure()
    # Plot the data
    map.plot(cmap=cmap, norm=norm)
    plt.title(var)
    # Show the plot
    plt.savefig(f'_figures/{var}_map.png')

# Perform pairwise diff between each possible combo of map
diff_pairs = list(product(maps, repeat=2))
diff_pairs = [(a, b) for a, b in diff_pairs if a != b]

for map1, map2 in diff_pairs:
    diff = map2.values - map1.values
    diff = xr.DataArray(diff, coords=map1.coords, dims=map1.dims)
    plt.figure()
    # Plot the data
    diff.plot(cmap='bwr', vmin=-500, vmax=500)
    plt.title(f'{map2.name} - {map1.name}')
    # Show the plot
    plt.savefig(f'_figures/{map2.name}_vs_{map1.name}_diff_map.png')



