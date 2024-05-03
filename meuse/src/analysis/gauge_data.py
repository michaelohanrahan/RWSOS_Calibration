from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from pandas import Timestamp

line_kw = {
	"color":"",
	"label":"",
	"linewidth":1.0
	}


def plot_gauge(
    xdata,
    ydata,
    title,
):
    """_summary_"""
    fig = plt.figure()
    fig.set_size_inches(8,5)
    ax = fig.add_subplot(111)
    ax.set_position([0.07, 0.07, 0.91, 0.85])

    kw = line_kw.copy()
    kw.update({"color": "#0000FF", "label": "Observed discharge"})
    ax.plot(xdata, ydata, **kw)

    ax.set_xlim([xdata[0], xdata[-1]])
    ax.set_ylim(0)
    ax.tick_params(axis="both", which="both", **{"labelsize": 7})
    ax.set_xlabel("Date", **{"fontsize":7})
    ax.set_ylabel("Discharge [m3/ s]", **{"fontsize":7})
    ax.set_title(title, **{"fontsize": 10})
    
    return fig


def plot_hist(
    hist: dict,
    title: str,
):
    """_summary_"""

    fig = plt.figure()
    fig.set_size_inches(8,5)
    ax = fig.add_subplot(111)
    ax.set_position([0.07, 0.07, 0.91, 0.85])

    ax.bar(hist.keys(), hist.values(), width=0.75, color='b')
    ax.tick_params(axis="both", which="both", **{"labelsize": 7})
    ax.set_xlabel("Year", **{"fontsize":7})
    ax.set_ylabel("Quantity", **{"fontsize":7})
    ax.set_title(title, **{"fontsize": 10})

    return fig


def deter_usable(n, q):
    if q == "Bad":
        return "Bad"
    if n == 0:
        return "Bad"
    if (q == "Medium" or q == "Weird") and n == 1:
        return "Bad"
    if (q == "Medium" and n >= 2):
        return "Medium"
    if (q == "Good" and n >= 2):
        return "Good"
    if (q == "Weird" and n >= 2):
        return "Maybe"


def sort_gauges(
    ds: xr.Dataset,
    site: gpd.GeoDataFrame,
    model: str,
    period: tuple,
    spinup: int,
    out: Path | str,
):
    """_summary_"""

    site["nyear"] = np.nan

    hist = {}

    for var in ds.data_vars:
        site_num = int(var)
        da = ds[var][ds[var].notnull()]
        sy = Timestamp(da.time[0].values).year
        ey = Timestamp(da.time[-1].values).year

        nyear = max(
            min(ey, period[1]) - max(sy, period[0]),
            0,
        )

        nyear = max(nyear - spinup, 0)        

        if site_num in site.SiteNumber.values:
            fig = plot_gauge(
                da.time.values,
                da.values,
                f"{da.site_name} ({da.name})"
            )

            for _y in range(sy, ey+1):
                if _y not in hist:
                    hist[_y] = 1
                    continue 
                hist[_y] += 1

            idx = site[site["SiteNumber"]==site_num].index[0]
            site.loc[idx, "nyear"] = nyear

            fig.savefig(
                Path(out, "visual", model, f"{var}.png"),
                dpi=300
            )
            plt.close()

            pass
        else:
            print(f"Gauge with number: {var} has not modelled data")
    
    site = site[site["nyear"].notna()]

    fig_hist = plot_hist(
        hist,
        model,
    )

    fig_hist.savefig(
        Path(out, "visual", f"{model}_hist.png"),
        dpi=300
    )

    site["Usability"] = site.apply(
        lambda row: deter_usable(row.nyear, row.Quality), axis=1
    )

    site.to_file(
        Path(out, f"{model}_obs.gpkg"),
    )


if __name__ == "__main__":
    ds = xr.open_dataset(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\discharge_obs.nc")
    kg = gpd.read_file(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\king_quality.gpkg")
    pg = gpd.read_file(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\pierce_quality.gpkg")
    out = Path(r"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES")
    sort_gauges(
        ds,
        kg,
        "king",
        (2010, 2019),
        1,
        out,
    )
    pass

