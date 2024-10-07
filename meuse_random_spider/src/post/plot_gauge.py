import os
import random
import sys
from datetime import datetime
from pathlib import Path

import matplotlib.dates as dts
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr

line_kw = {
	"color":"",
	"label":"",
	"linewidth":1.0
	}


def read_gauge(
    p: Path | str,
):
    """_summary_."""
    x = pd.read_csv(p, parse_dates=True, index_col="datetime")
    x.index -= pd.DateOffset(minutes=15)
    x = x["Q"].resample("H").mean()
    x = x * (0.3048**3)
    return x


def read_model_data(
    p: Path | str,
    gauge: str,
):
    """_summary_."""
    z = xr.open_dataset(p)
    vals = z.Q.sel(time=slice("2016-01-01","2017-12-31",None), Q_gauges_obs=gauge)
    return vals


def plot(
    gauge: str,
    results: tuple | list,
    obs: pd.DataFrame,
):
    """_summary_."""
    series = ["24hourly", "12hourly", "6hourly", "3hourly"]
    results = [
        item for _s in series for item in results if str(item).endswith(_s)
    ] + [results[-1]]
    colors = ["#FF0000", "#FFA200", "#04FF00", "#002BFF", "#00F0FF", "#007676"]
    # colors += [
    #     "#" + "".join([f"{random.randint(0,255):02x}" for _ in range(3)]) for _ in 
    #     range(len(results) - 1)
    # ]
    fig = plt.figure()
    fig.set_size_inches(w=16, h=10)
    ax = fig.add_subplot(111)
    ax.set_position([0.07, 0.07, 0.91, 0.85])

    for idx, res in enumerate(results):
        ds = read_model_data(Path(res, "output_scalar.nc"), gauge)
        md_kw = line_kw.copy()
        md_kw.update({"color": colors[idx], "label": f"Modelled discharge {res.stem}"})
        ax.plot(ds.time.values, ds.values, **md_kw)
    obs_kw = line_kw.copy()
    obs_kw.update({"color": "#000000", "label": "Observed discharge"})
    dates = obs.index
    ax.plot(dates, obs.values, **obs_kw)

    ax.set_xlim([datetime(2016,5,1), datetime(2016,7,15)])
    ax.set_ylim([0,200])
    ax.tick_params(axis="both", which="both", **{"labelsize": 7})
    ax.set_xlabel("Date", **{"fontsize":7})
    ax.set_ylabel("Discharge [m3/ s]", **{"fontsize":7})
    ax.set_title(f"Discharge at gauge: {gauge}", **{"fontsize": 10})
    ax.grid(linestyle="--")

    ax.legend(bbox_to_anchor=[0.01,0.98],loc='upper left', fontsize=8)

    pass


if __name__ == "__main__":
    obs = r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\12101500_15m_2Y.csv"
    obs = read_gauge(obs)
    out_dir = Path(r"p:\1000365-002-wflow\tmp\usgs_wflow\models\MODELDATA_PIERCE_500M\run_default")
    results = [
        item for item in out_dir.iterdir()
    ]
    gauge = "12101500"

    plot(
        gauge,
        results,
        obs
    )
    pass