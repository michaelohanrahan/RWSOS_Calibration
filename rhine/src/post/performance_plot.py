import xarray as xr
import os 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
import json

os.chdir('p:/11209265-grade2023/wflow/RWSOS_Calibration')
basin = 'meuse'
param_dims = ['f', 'ksat', 'ml', 'nr', 'rd', 'st']
levels = ['level0', 'level1']
with open(f'{basin}/config/MINIMAL_calib_recipe.json') as f:
    recipe = json.load(f)
lname = recipe.keys()
sname = [recipe[l]['short_name'] for l in lname]
method_t = {'mult':'scale factor', 'add':'offset value'}
method = [method_t[recipe[l]['method']] for l in lname]
s_to_l = dict(zip(sname, lname))
s_to_m = dict(zip(sname, method))

font_title = {'family': 'serif',
              'color':  'black',
              'weight': 'normal',
              'size': 16,
              }
font_label = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 14,
                }

#%% Plotting performance metrics
def plot_param_vs_metric():
    for level in levels:
        out_dir = f'{basin}/data/2-interim/calib_data/{level}/visualisation/'
        file =  f'{basin}/data/2-interim/calib_data/{level}/performance.zarr/'
        ds = xr.open_zarr(file)
        print(F"Processing {level}")
        for var in ds.data_vars:
            for dim in param_dims:
                fig, ax = plt.subplots()
                for gauge in ds.gauges.values:
                    dg = ds[var].sel(gauges=gauge)
                    df = dg.to_dataframe()
                    df = df.reset_index()
                    s1 = df[var]
                    s2 = df[dim]
                    ax.scatter(s1, s2)
                if var in s_to_l:
                    ax.set_xlabel(f"{s_to_l[var].capitalize()} {s_to_m[var]}", fontdict=font_label)
                else:
                    ax.set_xlabel(var.capitalize(), fontdict=font_label)
                if '_' in dim:
                    ax.set_ylabel(' '.join(dim.split('_')).capitalize(), fontdict=font_label)
                else:
                    ax.set_ylabel(dim.capitalize(), fontdict=font_label)
                ax.set_title(f'{var} vs {dim}, level: {level[-1]}'.capitalize(), fontdict=font_title)
                plt.savefig(f'{out_dir}/{var}_vs_{dim}.png')
                plt.close()

plot_param_vs_metric()            

#%% Plotting error surface
def plot_error_surface():
    for level in levels:
        out_dir = f'{basin}/data/2-interim/calib_data/{level}/visualisation/'
        file = f'{basin}/data/2-interim/calib_data/{level}/performance.zarr/'
        ds = xr.open_zarr(file)
        print(f"Processing {level}")
        param_pairs = list(itertools.combinations(param_dims, 2))

        # Plot each pair of parameters vs. the error values
        for (param1, param2) in param_pairs:
            for var in ds.data_vars:
                for gauge in ds.gauges.values:
                    dg = ds[var].sel(gauges=gauge)
                    df = dg.to_dataframe().reset_index()
                    
                    x = df[param1]
                    y = df[param2]
                    z = df[var]
                    
                    fig, ax = plt.subplots()
                    heatmap = ax.tricontourf(x, y, z, levels=14, cmap='viridis')
                    fig.colorbar(heatmap, ax=ax, label=var.capitalize())
                    
                    if '_' in param1:
                        ax.set_xlabel(' '.join(param1.split('_')).capitalize(), fontdict=font_label)
                    else:
                        ax.set_xlabel(param1.capitalize(), fontdict=font_label)
                    
                    if '_' in param2:
                        ax.set_ylabel(' '.join(param2.split('_')).capitalize(), fontdict=font_label)
                    else:
                        ax.set_ylabel(param2.capitalize(), fontdict=font_label)
                    
                    ax.set_title(f'{var} vs {param1} and {param2}, gauge: {gauge}, level: {level[-1]}'.capitalize(), fontdict=font_title)
                    outplot = f'{out_dir}/{gauge}'
                    os.makedirs(outplot, exist_ok=True)
                    plt.tight_layout()
                    plt.savefig(os.path.join(outplot, f'{var}_vs_{param1}_and_{param2}.png'), dpi=300)
                    plt.close()

plot_error_surface()     