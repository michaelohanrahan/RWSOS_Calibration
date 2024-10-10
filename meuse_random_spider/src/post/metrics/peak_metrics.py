"""
Created on Jan 12, 2024

@author: Jing Deng

This python script develops metrics to evaluate peak timing and/or magnitude 
from the model results (sim) compared to measurement data (obs)

"""

import numpy as np
import pandas as pd
from scipy import stats, signal
import xarray as xr
from xarray.core.dataarray import DataArray
import matplotlib.pyplot as plt

import os
import sys


def peak_timing_errors(obs: DataArray,
                sim: DataArray,
                window: int, 
                distance: int = None,
                prominence: float = None, 
                datetime_coord: str = None):
    """Difference in peak flow timing.
    Uses scipy.find_peaks to find peaks in the observed time series. Starting with all observed peaks, those with a
    prominence of less than half of standard deviation of the observed time series are discarded. And the lowest peaks
    are subsequently discarded until all remaining peaks have a distance of at least 24*3 steps. Finally, the
    corresponding peaks in the simulated time series are searched in a window of size `window` on either side of the
    observed peaks and the absolute time differences between observed and simulated peaks is calculated.
    
    Parameters
    ----------
    obs : DataArray
        Observed time series.
    sim : DataArray
        Simulated time series.
    window : int
        Size of window to consider on each side of the observed peak for finding the simulated peak. That is, the total
        window length to find the peak in the simulations is :math:`2 * \\text{window} + 1` centered at the observed
        peak.
    distance: int, optional
        Required minimal horizontal distance (>= 1) in samples between neighbouring peaks. 
        Larger distance will filter out peaks that are close in time.
        Default value is 24*3=72
    prominence: float, optional
        Required prominence of peaks. The peaks with a prominence less than this are discarded.
        Larger prominence will filter out peaks that are close in time and magnitude.
        Default value is np.std(obs.values).
    datetime_coord : str, optional
        Name of datetime coordinate. Tried to infer automatically as 'time' if not specified.

    Returns
    -------
    peaks : numpy array (datetime64)
        Datetime indices of peaks in obs. 
    timing_errors : list (float)
        Difference in peak timing. Positive value indicates simulated peak is late. Negative value indicates simulated peak is early.

    References
    -------
    https://github.com/neuralhydrology/neuralhydrology/blob/master/neuralhydrology/evaluation/metrics.py#L538
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    
    """   
    
    if distance is None:
        distance = 24*3  # default value as 24*3
    
    if prominence is None:
        prominence = np.std(obs.values)  # default value as 0.5 * np.std(obs.values)
        
    if datetime_coord is None:
        datetime_coord = 'time'  # default value as 'time'
    
    # get indices of peaks and their corresponding height
    peaks,_ = signal.find_peaks(obs.values, distance=distance, prominence=prominence)
    
    sim = sim.set_index({datetime_coord:'time'})
    obs = obs.set_index({datetime_coord:'time'})
    
    # convert peak indices to datetime indices
    peaks = obs[datetime_coord].values[peaks]
    
    # evaluate timing
    valid_peaks = []
    timing_errors = []
    window = pd.Timedelta(hours=int(window))
    
    
    
    for idx in peaks:
        # To make peaks datetime stamps we need to make window a timedelta
        # skip peaks at the start and end of the sequence and peaks around missing observations
        
        if (idx - window < sim.time.min()) or (idx + window > sim.time.max()) or (pd.date_range(start=idx - window, end=idx + window, freq='1H').size != window/pd.Timedelta(hours=1)*2 + 1):
            continue
        valid_peaks.append(idx)

        # check if the value at idx is a peak (both neighbors must be smaller)
        if (sim.loc[idx] > sim.loc[idx - pd.Timedelta(hours=1)]) and (sim.loc[idx] > sim.loc[idx + pd.Timedelta(hours=1)]):
            peak_sim = sim.loc[idx]
        else:
            # define peak around idx as the max value inside of the window
            values = sim.loc[idx - window : idx + window]
            if not values.isnull().all():
                peak_sim = values[values.argmax()]
            else:
                # Handle the case when all values are NaN
                peak_sim = np.nan

        # If peak_sim is NaN, skip this iteration
        if pd.isnull(peak_sim):
            timing_errors.append(np.nan)
            continue

        # get xarray object of qobs peak, for getting the date and calculating the datetime offset
        peak_obs = obs.loc[idx]

        # calculate the time difference between the peaks (positive value: sim is late; negative value: sim is early)
        delta = peak_sim.time - peak_obs.time
        timing_error = delta.values / pd.to_timedelta('1H')
        timing_errors.append(timing_error)

    return np.array(valid_peaks), timing_errors


# def peak_timing_errors(obs: DataArray,
#                 sim: DataArray,
#                 window: int, 
#                 distance: int = None,
#                 prominence: float = None, 
#                 datetime_coord: str = None):
#     """Difference in peak flow timing.
#     Uses scipy.find_peaks to find peaks in the observed time series. Starting with all observed peaks, those with a
#     prominence of less than half of standard deviation of the observed time series are discarded. And the lowest peaks
#     are subsequently discarded until all remaining peaks have a distance of at least 24*3 steps. Finally, the
#     corresponding peaks in the simulated time series are searched in a window of size `window` on either side of the
#     observed peaks and the absolute time differences between observed and simulated peaks is calculated.
    
#     Parameters
#     ----------
#     obs : DataArray
#         Observed time series.
#     sim : DataArray
#         Simulated time series.
#     window : int
#         Size of window to consider on each side of the observed peak for finding the simulated peak. That is, the total
#         window length to find the peak in the simulations is :math:`2 * \\text{window} + 1` centered at the observed
#         peak.
#     distance: int, optional
#         Required minimal horizontal distance (>= 1) in samples between neighbouring peaks. 
#         Larger distance will filter out peaks that are close in time.
#         Default value is 24*3=72
#     prominence: float, optional
#         Required prominence of peaks. The peaks with a prominence less than this are discarded.
#         Larger prominence will filter out peaks that are close in time and magnitude.
#         Default value is np.std(obs.values).
#     datetime_coord : str, optional
#         Name of datetime coordinate. Tried to infer automatically as 'time' if not specified.

#     Returns
#     -------
#     peaks : numpy array (int)
#         idx of peaks in obs. 
#     timing_errors : list (float)
#         Difference in peak timing. Positive value indicates simulated peak is late. Negative value indicates simulated peak is early.

#     References
#     -------
#     https://github.com/neuralhydrology/neuralhydrology/blob/master/neuralhydrology/evaluation/metrics.py#L538
#     https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
    
#     """   
    
#     if distance is None:
#         distance = 24*3  # default value as 24*3
    
#     if prominence is None:
#         prominence = np.std(obs.values)  # default value as 0.5 * np.std(obs.values)
        
#     if datetime_coord is None:
#         datetime_coord = 'time'  # default value as 'time'
    
#     # get indices of peaks and their corresponding height
#     peaks,_ = signal.find_peaks(obs.values, distance=distance, prominence=prominence)
    
#     # evaluate timing
#     valid_peaks = []
#     timing_errors = []
#     for idx in peaks:
#         # skip peaks at the start and end of the sequence and peaks around missing observations
#         if (idx - window < 0) or (idx + window >= len(obs)) or (pd.date_range(obs[idx - window][datetime_coord].values,
#                                                                               obs[idx + window][datetime_coord].values,
#                                                                               freq='1H').size != 2 * window + 1):
            
#             continue
        
#         valid_peaks.append(idx)

#         # check if the value at idx is a peak (both neighbors must be smaller)
#         if (sim[idx] > sim[idx - 1]) and (sim[idx] > sim[idx + 1]):
#             peak_sim = sim[idx]
        
#         else:
#             # define peak around idx as the max value inside of the window
#             #TODO: experiment to see if we can select from within the window in another way
            
#             values = sim[idx - window:idx + window + 1]
            
#             if not values.isnull().all():
#                 peak_sim = values[values.argmax()]
                
#             else:
#                 # Handle the case when all values are NaN
#                 # For example, you could set peak_sim to NaN or to a specific value
#                 peak_sim = np.nan
        
#         # If peak_sim is NaN, skip this iteration
#         if pd.isnull(peak_sim):
#             timing_errors.append(np.nan)
#             continue

#         # get xarray object of qobs peak, for getting the date and calculating the datetime offset
#         peak_obs = obs[idx]

#         # calculate the time difference between the peaks (positive value: sim is late; negative value: sim is early)
#         delta = peak_sim.coords[datetime_coord] - peak_obs.coords[datetime_coord]
#         timing_error = delta.values / pd.to_timedelta('1H')
#         timing_errors.append(timing_error)
    
#     return np.array(valid_peaks), timing_errors




# if __name__ == '__main__':
    
#     # set working folder
#     working_folder=r'p:/11209265-grade2023/wflow/wflow_meuse_julia/wflow_meuse_202401'
#     sys.path.append(working_folder)
    
#     # load example data: dataset that contains both obs. and model runs results
#     fn_ds = r'/_output/ds_obs_model_combined.nc'
#     ds = xr.open_dataset(working_folder+fn_ds)
    
#     # get the obs for Chooz (wflow_id=4)
#     obs = ds.sel(runs='Obs.', wflow_id=4).Q
#     print(type(obs))
    
#     # get the model run s07 results for Chooz (wflow_id=4)
#     sim = ds.sel(runs='s07', wflow_id=4).Q
#     print(type(sim))
    
#     # compute peak_timing_errors
#     peaks, timing_errors = peak_timing_errors(obs, sim, window=72)
    
    
#     # some example post analysis of the peak timing errors results:
#     # 1) mean peak timing: mean absolute time difference across all peaks
#     mean_peak_timing = np.mean(np.abs(timing_errors)) if len(timing_errors) > 0 else np.nan
    
#     # 2) mean absolute percentage peak error
#     obs_Q = obs[peaks].values
#     sim_Q = sim[peaks].values
#     peak_mape = np.sum(np.abs((sim_Q - obs_Q) / obs_Q)) / peaks.size * 100
    
#     # 3) scatter plot of timing_errors vs. obs
#     plt.scatter(obs[peaks], timing_errors)
#     plt.xlabel('obs Q, m3/s')
#     plt.ylabel('timing error, h')