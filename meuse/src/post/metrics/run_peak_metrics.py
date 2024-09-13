from metrics.peak_metrics import peak_timing_errors
import numpy as np
import pandas as pd
from icecream import ic



def store_peak_info(ds, id_key, window):
    """
    Store peak timing information in a dictionary.

    Parameters:
    - ds (xarray.Dataset): Dataset containing the data.
    - df_GaugeToPlot (pandas.DataFrame): DataFrame containing the gauge information.
    - id_key (str): Key to identify the station in the dataset.
    - window (int): Window size for peak timing errors.

    Returns:
    - peak_dict (dict): Dictionary containing the peak timing information for each station and run.
    """
    peak_dict = {}

    for id in ds[id_key].values:
        station_id = id

        # select a station using the id grouping and the sub selection id (station_id)
        ds_sub = ds.sel({id_key:station_id})

        # get obs data
        obs = ds_sub.sel(runs='Obs.').Q
        obs = obs.ffill('time')

        peak_dict[id] = {}

        for run in ds_sub.runs.values:
            if run != 'Obs.':
                
                #Sometimes Q is empty, we'd rather see what is empty in the plot than have an error
                if ds_sub.sel(runs=run).Q.isnull().all():
                    continue
                
                sim = ds_sub.sel(runs=run).Q
                
                peaks, timing_errors = peak_timing_errors(obs, sim, window=window)

                # print(f'run, peaks, timing_errors: {run}, \n {peaks}, \n  {timing_errors}')
                # Check if peaks is empty
                if len(peaks) > 0 and not np.isnan(peaks).all() and not np.isnan(timing_errors).all():
                    
                    timing_errors = np.array(timing_errors)
                    peaks = pd.Series(peaks)
                
                    #some timing errors should be nan where they aree the same as the window.
                    mask = np.abs(timing_errors) == window
                    
                    if np.any(mask):
                        # print('MASK', mask)
                        # print('where', np.where(mask)[0])
                        timing_errors[np.where(mask)[0]] = np.nan
                        peaks[np.where(mask)[0]] = pd.NaT
                    
                    peaks = np.array(peaks)
                    
                    # Convert timing_errors to timedelta (assuming timing_errors are in hours)
                    timing_errors_timedelta = pd.to_timedelta(timing_errors, unit='h')
                    
                    # print('timing_errors', timing_errors)   
                    # print('peaks', peaks)
                    
                    # Add timedelta to datetime
                    peaks_sim = peaks + timing_errors_timedelta

                    mean_peak_timing = np.nanmean(np.abs(timing_errors))
                    
                    peaks_index = pd.DatetimeIndex(peaks).dropna()
                    
                    obs_Q = obs.sel(time=peaks_index).values
                    sim_Q = sim.sel(time=peaks_index).values
                    
                    peak_mape = np.sum(np.abs((sim_Q - obs_Q) / obs_Q)) / peaks.size * 100
                
                else:
                    peaks_sim = np.nan
                    mean_peak_timing = np.nan
                    peak_mape = np.nan

                # Expand the inner dictionary with the inner loop
                peak_dict[id][run] = {'peaks': peaks_sim.dropna(), 
                                      'timing_errors': timing_errors[~np.isnan(timing_errors)], 
                                      'mean_peak_timing': mean_peak_timing,
                                      'peak_mape': peak_mape}

        peaks_obs, _ = peak_timing_errors(obs, obs, window=window)  # Calculate peaks for 'Obs.' with itself
        peak_dict[id]['Obs.'] = {'peaks': peaks_obs, 
                                 'timing_errors': np.zeros_like(peaks_obs), 
                                 'mean_peak_timing': np.nan,
                                 'peak_mape': np.nan}

    return peak_dict