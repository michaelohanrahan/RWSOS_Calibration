import xarray as xr 
import pandas as pd
import os 
from argparse import ArgumentParser
import pathlib


ds = xr.open_dataset(r'discharge_obs_HBV_combined.nc')

if os.path.exists('wflow_ids_in_obs.txt'):
    os.remove('wflow_ids_in_obs.txt')

def loop_data(ds, id_str, data):
    for i in ds[id_str].values:
        with open(f'{id_str}_in_{data}.txt', 'a') as f:
            len_series = len(ds.sel(wflow_id=i).time)
            num_na = ds.sel(wflow_id=i).Q.isnull().sum().values
            num_zeros = (ds.sel(wflow_id=i).Q == 0).sum().values
            freq = pd.infer_freq(ds.sel(wflow_id=i).time.values)
            pct_na = num_na/len_series
            f.write('ID: '+ str(i) + ' length: ' + str(len_series) + ' NA: ' + str(num_na) + ' Zeros: ' + str(num_zeros) + ' Freq: ' + str(freq) + '\n')


if "__name__" == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--id_str', type=str, default='wflow_id')
    parser.add_argument('--data', type=str, default='Q')
    parser.add_argument('--path', type=str, default=os.getcwd()) # default to
    
    id_str = parser.parse_args().id_str
    data_str = parser.parse_args().data
    cwd = parser.parse_args().path.as_posix()
    os.chdir(cwd)
    loop_data(ds, id_str, data_str)