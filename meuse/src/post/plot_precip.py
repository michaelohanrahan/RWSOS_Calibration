"""Plotting of gridded precipitation."""
from itertools import product
from pathlib import Path

import matplotlib.dates as dns
import matplotlib.lines as lns
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from hydromt import DataCatalog

from puget.plot import set_ax_props
from puget.util import CMIP6_MODELS, check_directory

SUBCATCH = {
    "KING": ["11", "21"],
    "PIERCE": ["20", "28", "40"],
}

SUBCATCH_RIVER = {
    "11": "Cedar River",
    "21": "Green River",
    "20": "Puyallup River",
    "28": "Chambers Creek",
    "40": "Nisqually River",
}


def yearly_max(
    da: xr.DataArray,
    sel: tuple | list,
    sel_id: str,
    periods: int,
):
    """_summary_."""
    da = da.sel({sel_id: sel})
    da = da.rolling(time=periods, center=True).sum()
    da = da.groupby("time.year").max()
    return da


def select_data(
    p: Path | str,
    tar_path: Path | str,
    corr_anno: Path | str,
    subcatch: tuple | list,
    periods: int,
):
    """_summary_."""
    tar = xr.open_dataset(
        Path(p, tar_path, "output_scalar.nc")
    )
    tar_corr = xr.open_dataset(
        Path(p, f"{tar_path}_{corr_anno}", "output_scalar.nc")
    )

    tar = yearly_max(tar.P, subcatch, "P_subcatchment", periods)
    tar_corr = yearly_max(tar_corr.P, subcatch, "P_subcatchment", periods)
    
    return tar, tar_corr


def plot_cdf(
    cl_name: str,
    river_name: str,
    tar_hist: xr.DataArray,
    tar_hist_corr: xr.DataArray,
    tar_fut: xr.DataArray,
    tar_fut_corr: xr.DataArray,
    output_dir: Path | str,
):
    """_summary_."""
    fig = plt.figure()
    fig.set_size_inches(8.3, 5.8)
    ax = fig.add_subplot(111)
    ax.set_position([0.065, 0.07, 0.92, 0.88])

    ax = set_ax_props(
        ax, 
        {"labelsize": 7}, 
        title=f"{river_name}; {cl_name}",
        title_props={"fontsize": 11.5}
    )

    maxx = 1

    p = (np.arange(tar_hist.year.size) + 1) / tar_hist.year.size
    rp = np.sort(1 / p)
    tar_hist_values = tar_hist.values
    tar_hist_corr_values = tar_hist_corr.values
    miny = min(np.min(tar_hist_values), min(tar_hist_corr_values))
    maxy = max(np.max(tar_hist_values), np.max(tar_hist_corr_values))
    maxx = max(maxx, len(rp))
    ax.plot(rp, np.sort(tar_hist_values), color="b", label="Historic")
    ax.plot(rp, np.sort(tar_hist_corr_values), color="b", linestyle="--", label="Historic (bias corrected)")

    p = (np.arange(tar_fut.year.size) + 1) / tar_fut.year.size
    rp = np.sort(1 / p)
    tar_fut_values = tar_fut.values
    tar_fut_corr_values = tar_fut_corr.values
    miny = min(miny, min(np.min(tar_fut_values), np.min(tar_fut_corr_values)))
    maxy = max(maxy, max(np.max(tar_fut_values), np.max(tar_fut_corr_values)))
    maxx = max(maxx, len(rp))
    ax.plot(rp, np.sort(tar_fut_values), color="r", label="Future")
    ax.plot(rp, np.sort(tar_fut_corr_values), color="r", linestyle="--", label="Future (bias corrected)")   

    ax.set_xlim(1,maxx)
    ax.set_ylim(miny,maxy*1.2)

    ax.set_xlabel("Return period (emprical) [Y]", fontsize=8)
    ax.set_ylabel("Precipitation (3 day avg.) [mm]", fontsize=8)

    ax.legend(bbox_to_anchor=(0.005, 0.99), loc="upper left", fontsize=8)

    fig.savefig(
        Path(output_dir, f"{cl_name}_{river_name.replace(' ', '_').lower()}.png"),
        dpi=300,
    )
    pass


if __name__ == "__main__":
    for region, cl_name in product(["KING", "PIERCE"], CMIP6_MODELS): 
        model_dir = Path(
            f"p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_{region}_CLIMATE"
        )
        tar_hist, tar_hist_corr = select_data(
            Path(model_dir, "run_default"),
            f"{cl_name}_historic",
            "bc",
            SUBCATCH[region],
            periods=24,
        )
        tar_fut, tar_fut_corr = select_data(
            Path(model_dir, "run_default"),
            f"{cl_name}_future",
            "bc",
            SUBCATCH[region],
            periods=24,
        )

        for sub in SUBCATCH[region]:
            river_name = SUBCATCH_RIVER[sub]
            plot_cdf(
                cl_name,
                river_name,
                tar_hist=tar_hist.sel(P_subcatchment=sub),
                tar_hist_corr=tar_hist_corr.sel(P_subcatchment=sub),
                tar_fut=tar_fut.sel(P_subcatchment=sub),
                tar_fut_corr=tar_fut_corr.sel(P_subcatchment=sub),
                output_dir="p:/1000365-002-wflow/tmp/usgs_wflow/docs/images/bias/cdf",
            )
    pass
