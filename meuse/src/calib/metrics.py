import hydromt.stats as stats
import numpy as np
import xarray as xr
from hydromt.stats import kge as kge_ds
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
    """_summary_"""
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
    """_summary_"""
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
    """_summary_"""
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
    """_summary_"""
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
    import xarray as xr
    obs = xr.open_dataset(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\discharge_obs_combined.nc")
    obs.Q.values = obs.Q.values * (0.3048**3)
    md = xr.open_dataset(r"p:\1000365-002-wflow\tmp\usgs_wflow\models\MODELDATA_KING_CALIB\calib_data\level1\ksat100_rd0.5_st0.5\run_default\output_scalar.nc")
    res = kge(
        md.sel(time=slice("2011-01-01", "2011-12-31")),
        obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
        gauges=["12119000", "12120000", "12113000"],
    )
    res2 = rld(
        md.sel(time=slice("2011-01-01", "2011-12-31")),
        obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
        gauges=["12119000", "12120000", "12113000"],
    )
    res3 = peakdis(
        md.sel(time=slice("2011-01-01", "2011-12-31")),
        obs.sel(time=slice("2011-01-01", "2011-12-31T00:00:00")),
        gauges=["12119000", "12120000", "12113000"],
    )
    weighted_euclidean(
        (np.array(res["kge"]), np.array(res2)),
        weights=(0.6, 0.4),
    )
    pass