from pathlib import Path
import numpy as np
import xarray as xr
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import pearsonr
import seaborn as sns
import plotly.graph_objects as go
import matplotlib.pyplot as plt

import dask
from dask import delayed


def smooth_observation_data(root_dir: str, 
                           observations_path: str, 
                           selected_gauge_path: str,
                           fig_dir: str,
                           plot_fig: list | None,
                           time_range: tuple = (None, None),
                           sigma: int = 5,
                           ):
    
    root = Path(root_dir)
    observations = root / observations_path
    selected_gauge = root / selected_gauge_path
    fig_path = Path(fig_dir)
    
    obs = xr.open_dataset(observations)
    df = pd.read_csv(selected_gauge)
    _temp = df.values.flatten()
    gauges = [int(num) for num in _temp if pd.notna(num)]
    
    # Extract the Q data for the selected gauges
    q_data = obs.sel(wflow_id=gauges, time=slice(*time_range))

    @delayed
    def process_gauge(gauge):
        da = q_data.sel(wflow_id=gauge).Q
        
        # interpolate/extrapolate the data to fill nan values maximun 3 hours
        da = da.interpolate_na(dim='time', method='linear', limit=3)
        
        _da_smooth = gaussian_filter1d(da, sigma=sigma, mode='nearest')
        
        da_smooth = xr.DataArray(
            np.float64(_da_smooth),
            dims=da.dims,
            coords=da.coords,
            name='Q',
            attrs=da.attrs
        )
        
        # plot to check the difference between da and da_smooth
        if plot_fig is not None and gauge in plot_fig:
            # get gauge name
            try:
                gauge_name = q_data.sel(wflow_id=gauge).station_names.item().decode('utf-8', errors='ignore')
            except:
                gauge_name = q_data.sel(wflow_id=gauge).station_names.item()
            # hydrograph check
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=da.time, y=da, mode='lines', name='Original Data'))
            fig.add_trace(go.Scatter(x=da_smooth.time, y=da_smooth, mode='lines', name='Smoothed Data'))
            fig.update_layout(title=f'{gauge_name}_{gauge}', xaxis_title='Time', yaxis_title='Discharge (m3/s)')
            fig.write_html(fig_path / f'hydrograph_compare/{gauge_name}_{gauge}.html')
            
            # Statistical measures
            mask = np.isfinite(da.values) & np.isfinite(da_smooth.values)
            da_clean = da.values[mask]
            da_smooth_clean = da_smooth.values[mask]
            correlation = pearsonr(da_clean, da_smooth_clean)[0]
            rmse = np.sqrt(np.mean((da_clean - da_smooth_clean)**2))

            # Plotting the density distribution of da and da_smooth
            fig, ax = plt.subplots()
            sns.kdeplot(da, ax=ax, label='Original Data')
            sns.kdeplot(da_smooth, ax=ax, label='Smoothed Data')
            ax.set_title(f'Density Distribution of Original vs Smoothed Data for {gauge_name}_{gauge}\nCorrelation: {correlation:.4f}, RMSE: {rmse:.4f}')
            ax.set_xlabel('Discharge (m3/s)')
            ax.set_ylabel('Density')
            plt.legend()
            plt.savefig(fig_path / f'density_distribution/{gauge_name}_{gauge}.png')
            plt.close()
            
        else:
            pass
        
        return da_smooth
        
    # Use Dask to process gauges in parallel
    delayed_results = [process_gauge(gauge) for gauge in gauges]
    das = dask.compute(*delayed_results)

    ds_smooth = xr.Dataset({'Q': xr.concat(das, dim='wflow_id')})

    # Initialize other vars in ds_smooth before assigning their values from q_data
    for var in ['x', 'y', 'z', 'lon', 'lat', 'station_id', 'station_names']:
        ds_smooth[var] = q_data[var]
    
    return ds_smooth



if __name__ == '__main__':
    
    root_dir = r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine'
    observations_path = 'data/1-external/discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727_with_new_data.nc'
    selected_gauge_path = 'data/2-interim/manually_selected_gauges_v2.csv'
    fig_dir = r'p:\11209265-grade2023\wflow\RWSOS_Calibration\rhine\data\5-visualization\smooth_obs'
    plot_fig = [688]
    time_range=('1996', None)
    sigma = 5
    
    import time
    start_time = time.time()
    ds_smooth = smooth_observation_data(root_dir, 
                                        observations_path, 
                                        selected_gauge_path, 
                                        fig_dir, 
                                        plot_fig, 
                                        time_range, 
                                        sigma)
    print(f"Time taken: {time.time() - start_time} seconds")
    
    ds_smooth.to_netcdf(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external\discharge_hourlyobs_smoothed.nc')
    
    
    
    # # check if obs vars are correctly copied to ds vars
    # obs = xr.open_dataset(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external\discharge_obs_hr_FORMAT_allvars_wflowid_0_to_727_with_new_data.nc')
    # time_range=('1996', None)
    # obs = obs.sel(time=slice(*time_range))
    
    # ds_smooth = xr.open_dataset(r'c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\rhine\data\1-external\discharge_hourlyobs_smoothed.nc')
    
    # obs_lob = obs.sel(wflow_id=709).Q
    # smooth_lob = ds_smooth.sel(wflow_id=709).Q
    
    # plt.figure(figsize=(10, 6))
    # # plt.plot(obs_lob.time, obs_lob, label='Original Data')
    # # plt.plot(smooth_lob.time, smooth_lob, label='Smoothed Data')
    # # plot the difference between obs and smooth
    # diff = obs_lob - smooth_lob
    # plt.plot(smooth_lob.time, diff, label='Difference')
    # plt.title('Lobith Discharge Comparison')
    # plt.xlabel('Time')
    # plt.ylabel('Discharge (m3/s)')
    # plt.legend()
    # plt.grid(True)
    # plt.show()
