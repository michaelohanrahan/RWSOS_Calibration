# Modified by Jing on July 17, 2024
# TO-DOs: 
## 1. do we need func: _peakdis, peakdis, fix_minmax, fix_maxmin, fix_gap?


import hydromt.stats as stats
from hydromt.stats import kge as kge_ds, skills
import numpy as np
import pandas as pd
import xarray as xr
from scipy.signal import argrelmax, argrelmin, find_peaks

from plot_rld import main as rld_fig


def kge(
    sim: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,
):
    """_summary_"""
    res = {
        "kge": [],
        "kge_pearson_coef": [],
        "kge_rel_var": [],
        "kge_bias": [],
    }

    for g in gauges:
        da = kge_ds(
            sim.sel(wflow_id=g).Q,
            obs.sel(wflow_id=g).Q,
        )

        for var in da.data_vars:
            res[var].append(
                round(float(da[var].values),4)
            )

    return res


def nse(
    sim: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,
):
    res = []
    
    for g in gauges:
        da = skills.nashsutcliffe(
            sim.sel(wflow_id=g).Q,
            obs.sel(wflow_id=g).Q,
        )
        res.append(round(float(da.values),4))
    
    return res


def nse_log(
    sim: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,
):
    res = []
    
    for g in gauges:
        da = skills.lognashsutcliffe(
            sim.sel(wflow_id=g).Q,
            obs.sel(wflow_id=g).Q,
        )
        res.append(round(float(da.values),4))
    
    return res
    

def mm7q(
    ds: xr.Dataset | xr.DataArray,
    dry_month: list,
):
    """
    Monthly minimum 7-day discharge for selected dry months.
    
    Args:
        ds (xr.Dataset | xr.DataArray): The dataset that contains the discharge time series for single station.
        dry_month (list): List of dry months. 
                          You can specify only the start and end months of dry period, e.g. [6,9]
                          Or you can specify each month, e.g., [6,7,8,9]
        
    Returns:
        xr.Dataset | xr.DataArray: Monthly minimum 7-day discharge for selected dry month.
    """
    # get the rolling mean window based on timescale in dataset
    if xr.infer_freq(ds.time).lower() == "D":
        window = 7
    elif xr.infer_freq(ds.time).lower() == "h":
        window = 7 * 24
    elif xr.infer_freq(ds.time).lower() == "3h":
        window = 7 * int(24 / 3)
    
    # calculate the MM7Q for each month
    _mm7q = ds.rolling(time=window).mean().resample(time='M').min('time').compute()
    
    # select out the MM7Q for selected dry months
    dry_month_start = dry_month[0]
    dry_month_end = dry_month[-1]
    months = _mm7q['time'].dt.month
    mm7q_dry_month = _mm7q.sel(time=_mm7q['time'].where((months>=dry_month_start)
                                                        &(months<=dry_month_end),
                                                        drop=True))
    
    return mm7q_dry_month
    
    
def nselog_mm7q(
    sim: xr.Dataset,
    obs: xr.Dataset,
    dry_month: list,
    gauges: tuple | list,
):
    """nse-log of mm7q of modeled discharge compared to observations for selected dry months and gauges

    Args:
        sim (xr.Dataset): Model dataset containing discharge values.
        obs (xr.Dataset): Observed dataset containing discharge values.
        dry_month (list): List of dry months.
        gauges (tuple | list): Tuple or list of gauges wflow_id for which needs to be calculated.

    Returns:
        List: List of nselog_mm7q values for each gauge (wflow_id).
    """
    
    res = []
    
    for g in gauges:
        sim_mm7q = mm7q(sim.sel(wflow_id=g).Q, dry_month)
        obs_mm7q = mm7q(obs.sel(wflow_id=g).Q, dry_month)
        nselog_mm7q = skills.lognashsutcliffe(sim_mm7q, obs_mm7q)
        res.append(round(float(nselog_mm7q.values),4))
    
    return res


def peak_timing_errors(
    sim: xr.DataArray,
    obs: xr.DataArray,
    window: int, 
    distance: int = None,
    prominence: float = None, 
    datetime_coord: str = None,
):
    """Difference in peak flow timing.
    Uses scipy.find_peaks to find peaks in the observed time series. Starting with all observed peaks, those with a
    prominence of less than half of standard deviation of the observed time series are discarded. And the lowest peaks
    are subsequently discarded until all remaining peaks have a distance of at least 24*3 steps. Finally, the
    corresponding peaks in the simulated time series are searched in a window of size `window` on either side of the
    observed peaks and the absolute time differences between observed and simulated peaks is calculated.
    
    Parameters
    ----------
    sim : xr.DataArray
        Simulated time series.
    obs : xr.DataArray
        Observed time series.
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
        prominence = np.nanstd(obs.values)  # default value as 0.5 * np.std(obs.values)
        
    if datetime_coord is None:
        datetime_coord = 'time'  # default value as 'time'
    
    # get indices of peaks and their corresponding height
    peaks,_ = find_peaks(obs.values, distance=distance, prominence=prominence)
    
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


def peak_errors(
    sim: xr.Dataset,
    obs: xr.Dataset,
    window: int,
    gauges: tuple | list,
):
    """peak error dictionary: mae of timing errors, and mape of peaks for gauges

    Args:
        sim (xr.Dataset): Model dataset containing discharge values.
        obs (xr.Dataset): Observed dataset containing discharge values.
        window (int): Size of window to consider on each side of the observed peak for 
                      finding the simulated peak.
        gauges (tuple | list): Tuple or list of gauges wflow_id for which needs to be calculated.

    Returns:
        Dict: Dictionary of mae_timing and mape_peak for each gauge (wflow_id).
    """
    
    res = {
        "mae_timing": [],
        "mape_peak": [],
    }
    
    for g in gauges:
        sim_g = sim.sel(wflow_id=g).Q
        obs_g = obs.sel(wflow_id=g).Q
        peaks , timing_errors = peak_timing_errors(sim_g, obs_g, window)
        
        # compute mae of timing_erros
        mae_timing = np.mean(np.abs(timing_errors))
        
        # compute mape of peak magnitude
        peaks_index = pd.DatetimeIndex(peaks).dropna()
        obs_peak = obs_g.sel(time=peaks_index).values
        sim_peak = sim_g.sel(time=peaks_index).values            
        peak_mape = np.sum(np.abs((sim_peak - obs_peak) / obs_peak)) / peaks.size * 100
                
        res['mae_timing'].append(round(float(mae_timing),4))
        res['mape_peak'].append(round(float(peak_mape),4))
    
    return res
    
    

def _rld(
    data: np.ndarray,
    thres: float = 95,
    buildup: bool = False,
    flow_based: bool = False,
):

    """
    Calculate the Rising Limb Density (RLD) metric.

    Args:
        data (np.ndarray): The input data array.
        thres (float, optional): The threshold percentile for filtering smaller events. Defaults to 95.
        buildup (bool, optional): Whether to group pre-events with larger events. Defaults to False.
        flow_based (bool, optional): Whether to calculate RLD based on flow. Defaults to False.

    Returns:
        float: The calculated RLD value.

    Raises:
        ValueError: If the input data is empty.

    """
    
    minx = find_peaks(data * -1)[0]
    maxx = find_peaks(data)[0]

    if len(maxx) == 0 or len(minx) == 0:
        return np.nan

    if maxx[0] < minx[0]: 
        maxx = maxx[1:]
    if len(maxx) == 0:
        return np.nan
    
    if maxx[-1] < minx[-1]:
        minx = minx[:-1]
    if len(minx) == 0:
        return np.nan
    
    # TODO look at this later on
    # Need to figure out how bad data can be 
    # and the influence on peak finding
    shape = (len(minx), len(maxx))
    shape_old = (np.nan,np.nan)
    while True:
        if shape == shape_old and shape[0] == shape[1]:
            break
        shape_old = shape
        if shape[0] > shape[1]:
            minx = fix_minmax(minx, maxx, shape)
        elif shape[1] > shape[0]:
            maxx = fix_maxmin(maxx, minx, shape)
        else:
            minx = fix_minmax(minx, maxx, shape)
            maxx = fix_maxmin(maxx, minx, shape)

        shape = (len(minx), len(maxx))
        if 0 in shape:
            return np.nan

    # Filter all smaller events
    temp = data[maxx] - data[minx]
    red = np.where(temp < np.percentile(temp, thres))[0]
    maxx = np.delete(maxx, red)
    minx = np.delete(minx, red)

    # Group pre-events with the larger events
    if buildup:
        idx = 0
        while idx < len(maxx):
            if idx == 0:
                idx += 1
                continue
            betw = data[maxx[idx - 1]] - data[minx[idx]]
            prev = data[maxx[idx - 1]] - data[minx[idx - 1]]

            q_prev = data[maxx[idx-1]] - data[minx[idx-1]]
            q_next = data[maxx[idx]] - data[minx[idx-1]]

            if (2 * betw < prev) and (q_next > 2.5 * q_prev):
                if sum(np.isnan(data[minx[idx-1]: minx[idx]])) > 0:
                    pass
                else:
                    maxx = np.delete(maxx, idx - 1)
                    minx = np.delete(minx, idx)
                    idx -= 1
            idx += 1

    res = maxx - minx

    if flow_based:
        q_diff = data[maxx] - data[minx]
        res = q_diff / (maxx - minx)
        return np.mean(res)
    
    return len(res) / sum(res)


def rld(
    sim: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,
):
    """_summary_"""
    res = []

    for g in gauges:
        sim_e = _rld(
            sim.sel(wflow_id=g).Q.values,
        )
        obs_e = _rld(
            obs.sel(wflow_id=g).Q.values,
        )
        if np.isnan(sim_e) or np.isnan(obs_e):
            e = 1
        else:
            e = 1 - abs(1 - (sim_e / obs_e))
        res.append(round(e, 4))

    return res


def _peakdis(
    data: tuple | list,
    upper: float | int = 90,
    lower: float | int = 50,
):
    """_summary_"""

    peaks = (
        np.nanpercentile(data, upper) - np.nanpercentile(data, lower)
    )
    e = peaks / (0.9 - 0.5)

    return e


def peakdis(
    sim: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,    
):
    """
    Calculate the peak discharge discrepancy between model (sim) and observed (obs) datasets for a given set of gauges.

    Parameters:
    sim (xr.Dataset): Model dataset containing discharge values.
    obs (xr.Dataset): Observed dataset containing discharge values.
    gauges (tuple | list): Tuple or list of gauge wflow_id for which the peak discharge discrepancy needs to be calculated.

    Returns:
    list: List of peak discharge discrepancies for each gauge.

    """
    res = []
    
    for g in gauges:
        sim_val = _peakdis(sim.sel(wflow_id=g).Q.values)
        obs_val = _peakdis(obs.sel(wflow_id=g).Q.values)

        if np.isnan(sim_val) or np.isnan(obs_val):
            r = 1
        else:
            r = 1 - abs(1 - (sim_val/ obs_val))
        res.append(r)

    return res


def fix_minmax(
    minx: tuple | list,
    maxx: tuple | list,
    shape: tuple | list,    
):
    """
    Adjusts the minx array to match the size of the maxx array based on the given shape.

    Args:
        minx (tuple | list): The minx array.
        maxx (tuple | list): The maxx array.
        shape (tuple | list): The desired shape of the minx array.

    Returns:
        tuple | list: The adjusted minx array.

    Raises:
        None

    """
    while True:
        size_diff = len(minx)-len(maxx)
        temp = maxx - minx[size_diff:]
        mask = temp < 0
        loc_rv = np.argmax(np.flip(mask))
        if loc_rv == 0:
            break
        
        loc = len(temp) - loc_rv - 1 + size_diff
        minx = np.delete(minx, loc)

        if shape[0] == shape[1]:
            return minx
        
        shape = (len(minx), len(maxx))
    return minx


def fix_maxmin(
    maxx: tuple | list,
    minx: tuple | list,
    shape: tuple | list,
):
    """
    Adjusts the `maxx` array by removing elements based on the `minx` array.
    
    Args:
        maxx (tuple | list): The maximum values array.
        minx (tuple | list): The minimum values array.
        shape (tuple | list): The shape of the arrays.
        
    Returns:
        tuple | list: The adjusted `maxx` array.
    """
    while True:
        temp = maxx[0:len(minx)] - minx
        loc = np.argmax(temp < 0)

        if loc == 0:
            break

        maxx = np.delete(maxx, loc)

        if shape[0] == shape[1]:
            return maxx
        
        shape = (len(minx), len(maxx))
    return maxx


def fix_gap():
    pass



def weighted_euclidean(
    coef: tuple | list,
    weights: tuple | list,
    weighted = True,
):
    """_summary_"""

    if weighted and len(weights) != len(coef) or sum(weights) != 1:
        raise ValueError("")
    
    if not weighted:
        weights = [1] * len(coef)

    dist = [
        w * (1-item)**2 for item, w in zip(coef, weights) 
    ]

    res = np.sqrt(sum(dist))

    return list(res.round(4))


if __name__ == "__main__":
    
    ds = xr.open_dataset(r'p:\11209265-grade2023\wflow\wflow_meuse_julia\wflow_meuse_20240529_flpN_landN\_output\ds_obs_model_combined.nc')
    sim = ds.sel(runs='scale_10')
    obs = ds.sel(runs='Obs.')
    
    kge_res = kge(sim, obs, [16, 801])
    nse_res = nse(sim, obs, [16, 801])
    nse_log_res = nse(sim, obs, [16, 801])
    
    nselog_mm7q_res = nselog_mm7q(sim, obs, [6,11], [16, 801])
    
    peak_res = peak_errors(sim, obs, 72, [16, 801])