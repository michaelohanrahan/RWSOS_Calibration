# Modified by Jing on July 1, 2024
# TO-DOs: 
# 1. change md to sim?
# 2. function nse_mm7q: make sure the gauge variable names (now using 'Q_gauges_obs') 
    # are aligned with the real output nc file 

import hydromt.stats as stats
from hydromt.stats import kge as kge_ds, skills
import numpy as np
import xarray as xr
from scipy.signal import argrelmax, argrelmin, find_peaks

from plot_rld import main as rld_fig


def kge(
    md: xr.Dataset,
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
            md.sel(Q_gauges_obs=g).Q,
            obs.sel(Q_gauges_obs=g).Q,
        )

        for var in da.data_vars:
            res[var].append(
                round(float(da[var].values),4)
            )

    return res


def nse():
    pass


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
    md: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,    
):
    """
    Calculate the peak discharge discrepancy between model (md) and observed (obs) datasets for a given set of gauges.

    Parameters:
    md (xr.Dataset): Model dataset containing discharge values.
    obs (xr.Dataset): Observed dataset containing discharge values.
    gauges (tuple | list): Tuple or list of gauge names for which the peak discharge discrepancy needs to be calculated.

    Returns:
    list: List of peak discharge discrepancies for each gauge.

    """
    res = []
    
    for g in gauges:
        md_val = _peakdis(md.sel(Q_gauges_obs=g).Q.values)
        obs_val = _peakdis(obs.sel(Q_gauges_obs=g).Q.values)

        if np.isnan(md_val) or np.isnan(obs_val):
            r = 1
        else:
            r = 1 - abs(1 - (md_val/ obs_val))
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
    md: xr.Dataset,
    obs: xr.Dataset,
    gauges: tuple | list,
):
    """_summary_"""
    res = []

    for g in gauges:
        md_e = _rld(
            md.sel(Q_gauges_obs=g).Q.values,
        )
        obs_e = _rld(
            obs.sel(Q_gauges_obs=g).Q.values,
        )
        if np.isnan(md_e) or np.isnan(obs_e):
            e = 1
        else:
            e = 1 - abs(1 - (md_e / obs_e))
        res.append(round(e, 4))

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
    
    
def nse_mm7q(
    md: xr.Dataset,
    obs: xr.Dataset,
    dry_month: list,
    gauges: tuple | list,
):
    """nse of mm7q of modeled discharge compared to observations for selected dry months and gauges

    Args:
        md (xr.Dataset): Model dataset containing discharge values.
        obs (xr.Dataset): Observed dataset containing discharge values.
        dry_month (list): List of dry months.
        gauges (tuple | list): Tuple or list of gauge names for which the peak discharge discrepancy needs to be calculated.

    Returns:
        List: List of nse_mm7q for each gauge.
    """
    # TO-DO: make sure the gauge variable names (now using 'Q_gauges_obs') 
    # are aligned with the real output nc file
    
    res = []
    
    for g in gauges:
        md_mm7q = mm7q(md.sel(Q_gauges_obs=g).Q, dry_month)
        obs_mm7q = mm7q(obs.sel(Q_gauges_obs=g).Q, dry_month)
        nse_mm7q = skills.nashsutcliffe(md_mm7q, obs_mm7q)
        res.append(nse_mm7q)
    
    return res
    


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
    # import xarray as xr
    # obs = xr.open_dataset(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\discharge_obs_combined.nc")
    # obs.Q.values = obs.Q.values * (0.3048**3)
     
    # res = kge(
    #     md.sel(time=slice("2011-01-01", "2011-12-31")),
    #     obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
    #     gauges=["12119000", "12120000", "12113000"],
    # )
    # res2 = rld(
    #     md.sel(time=slice("2011-01-01", "2011-12-31")),
    #     obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
    #     gauges=["12119000", "12120000", "12113000"],
    # )
    # res3 = peakdis(
    #     md.sel(time=slice("2011-01-01", "2011-12-31")),
    #     obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
    #     gauges=["12119000", "12120000", "12113000"],
    # )
    # weighted_euclidean(
    #     (np.array(res["kge"]), np.array(res2)),
    #     weights=(0.6, 0.4),
    # )
    
    # dumb test function mm7q
    obs = xr.open_dataset(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\1-external\discharge_obs_combined_EXAMPLE.nc')
    obs_sel = obs.sel(time=slice("2009-01-01", "2020-12-31"))
    
    res = nse_mm7q(obs_sel, obs_sel, [6,9], ['12105900'])