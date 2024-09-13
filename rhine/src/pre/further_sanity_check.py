from pathlib import Path
import dask.array as da
import numpy as np
import xarray as xr
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def monthly_nan_percentage(root_dir: str, 
                           observations_path: str, 
                           selected_gauge_path: str,
                           time_range: tuple = (None, None)):
    root = Path(root_dir)
    observations = root / observations_path
    selected_gauge = root / selected_gauge_path
    
    obs = xr.open_dataset(observations)
    df = pd.read_csv(selected_gauge)
    _temp = df.values.flatten()
    gauges = [int(num) for num in _temp if pd.notna(num)]
    
    # 1. Extract the Q data for the specified gauges
    q_data = obs.sel(wflow_id=gauges, runs='Obs.', time=slice(*time_range)).Q
    
    # 2. Compute the monthly nan percentage
    q_data_dask = q_data.chunk({'time': -1, 'wflow_id': 1})
    monthly_nan_percentage = q_data_dask.isnull().resample(time='1M').mean() * 100

    # 3. Convert to a pandas DataFrame for easier plotting
    df = monthly_nan_percentage.to_pandas().T
    return df


if __name__ == '__main__':
    root_dir = r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine'
    observations_path = 'data/1-external/discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727.nc'
    selected_gauge_path = 'data/2-interim/manually_selected_gauges.csv'
    
    df = monthly_nan_percentage(root_dir, observations_path, selected_gauge_path, time_range=('1996', None))

    # Create a heatmap
    plt.figure(figsize=(50, len(df.index)*0.3))  # Adjust figure height dynamically based on number of gauges
    sns.heatmap(df, cmap='YlOrRd', cbar_kws={'label': 'Missing Data Percentage', 'shrink': 0.5})

    plt.title('Monthly Missing Data Percentage for Selected Gauges')
    plt.xlabel('Year and Month')
    plt.ylabel('Gauge ID', rotation=0)

    plt.xticks(ticks=np.arange(len(df.columns)), labels=df.columns.strftime('%Y-%m'))
    
    # Add vertical lines at 2015-01 and 2016-01
    date_2015 = pd.to_datetime('2005-01-01')
    date_2016 = pd.to_datetime('2016-01-01')

    x_2015 = np.where(df.columns >= date_2015)[0][0]
    x_2016 = np.where(df.columns >= date_2016)[0][0]

    plt.axvline(x=x_2015, color='k', linestyle='--')
    plt.axvline(x=x_2016, color='k', linestyle='--')

    plt.show()
    
    
    # create hydrograph
    work_dir = Path(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine')
    output_dir = Path(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\rhine\data\5-visualization\observation_hydrograph')
    obs = xr.open_dataset(work_dir / observations_path)
    _df = pd.read_csv(work_dir / selected_gauge_path)
    _temp = _df.values.flatten()
    gauges = [int(num) for num in _temp if pd.notna(num)]
    
    time_range=('1996', None)
    q_data = obs.sel(wflow_id=gauges, runs='Obs.', time=slice(*time_range))
    
    # Plot the hydrograph
    from tqdm import tqdm
    for gauge in tqdm(gauges, desc="Processing gauges"):
        q_data_exp = q_data.sel(wflow_id=gauge)
        try:
            gauge_name = q_data_exp.station_names.item().decode('utf-8', errors='ignore')
        except:
            gauge_name = q_data_exp.station_names.item()
        
        plt.figure(figsize=(20, 6), dpi=300)
        plt.plot(q_data_exp.time, q_data_exp.Q, label=f'Gauge {gauge_name}_{gauge}')
        plt.title(f'Hydrograph for Gauge {gauge_name}_{gauge}')
        plt.xlabel('Time')
        plt.ylabel('Discharge (m3/s)')
        plt.legend()
        plt.savefig(output_dir / f'{gauge_name}_{gauge}.png', format='png')
        plt.close()

