import copy
import datetime
import math
import os
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import xclim
import yaml
from xclim import sdba
from xclim.sdba.adjustment import QuantileDeltaMapping
from hydromt import DataCatalog

from puget.io import export_to_zarr
from puget.util import (
    CMIP6_CATALOG_BP,
    CMIP6_MODELS,
    MyDumper,
    UNITS_MAP,
    UNITS_MULT, 
    check_directory, 
    make_path_relative,
)


def set_ax_props(
    ax: plt.axes,
    title: str = "bla",
):
    """_summary_."""
    ax.tick_params(axis="both", which="both", **{"labelsize": 7})
    ax.set_title(title, **{"fontsize": 9.5})
    return ax


def set_xyticks(
    ax: plt.axes,
    bounds: list | tuple,
):
    """_summary_."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xrange = xlim[1] - xlim[0]
    yrange = ylim[0] - ylim[1]
    xbounds = bounds[2] - bounds[0]
    ybounds = bounds[1] - bounds[3]

    xticks = list(
        range(math.ceil(bounds[0]), math.floor(bounds[2]) + 1)
    )
    yticks = list(
        range(math.ceil(bounds[1]), math.floor(bounds[3]) + 1)
    )

    xloc = [
        ((item - bounds[0])/xbounds)*xrange + xlim[0] for item in xticks
    ]
    yloc = [
        ((item - bounds[1])/ybounds)*yrange + ylim[0] for item in yticks
    ]

    ax.set_xticks(xloc)
    ax.set_yticks(yloc)
    ax.set_xticklabels(xticks)
    ax.set_yticklabels(yticks)
    return ax


def add_colorbar(
    fig: plt.figure,
    vrange: tuple,
    title: str,
    unit: str = "-",
):
    """_summary_"""
    cl_ax = fig.add_axes([0.75, 0.04, 0.02, 0.42])
    norm = mpl.colors.Normalize(vmin = vrange[0], vmax = vrange[1])
    cmp = mpl.colorbar.ColorbarBase(cl_ax, cmap="viridis", norm=norm)

    cl_ax.tick_params(axis="both", which="both", **{"labelsize": 8})
    cl_ax.set_ylabel(f"{title} [{unit}]", **{"fontsize": 9})
    return cl_ax, cmp 


def plot_correction(
    ref: np.ndarray,
    tar: np.ndarray,
    corr: np.ndarray,
    title: str,
):
    """_summary_."""
    fig = plt.figure()
    fig.set_size_inches(8,5)
    ax = fig.add_subplot(111)
    ax.set_position([0.05, 0.05, 0.93, 0.87])
    pass


def plot_spatial_correction(
    ref: xr.DataArray,
    tar: xr.DataArray,
    corr: xr.DataArray,
    title: str,
    scenario: str,
    quantile: float = 0.5,
):
    """_summary_."""
    # Figure
    fig = plt.figure()
    fig.set_size_inches(8,8)
    
    # Set the axes
    axes = []
    ax = fig.add_subplot(221)
    ax.set_position([0.05, 0.53, 0.42, 0.42])
    axes.append(ax)
    ax2 = fig.add_subplot(222)
    ax2.set_position([0.55, 0.53, 0.42, 0.42])
    axes.append(ax2)
    ax3 = fig.add_subplot(223)
    ax3.set_position([0.29, 0.04, 0.42, 0.42])
    axes.append(ax3)

    fig.text(x=0.5, y=0.99, s=f"{title} ({quantile} quantile)", ha="center", va="top", fontsize=12)
    
    oned = [
        item.quantile(quantile, dim="time").values for item in [ref, tar, corr]
    ]

    vmin = min([item.min() for item in oned])
    vmax = max([item.max() for item in oned])

    # Figure titles
    _titles = ["ERA5", f"{scenario}", f"{scenario}_corrected"]
    for idx, _ax in enumerate(axes):
        _ax.imshow(oned[idx], vmin=vmin, vmax=vmax, cmap="viridis", aspect="auto")
        _ax = set_ax_props(_ax, title=_titles[idx])
        _ax = set_xyticks(_ax, tar.raster.bounds)
    
    cl_ax, _ = add_colorbar(fig, [vmin, vmax], title, ref.units)
    return fig


def var_correction(
    tar: xr.DataArray,
    ds_map: xr.Dataset,
):
    """_summary_."""
    QDM = QuantileDeltaMapping.from_dataset(ds_map)
    da_corr = QDM.adjust(
        sim=tar,
        interp="linear",
    )
    return da_corr


def train_mapping(
    ref: xr.DataArray,
    tar: xr.DataArray,
    quantiles: list | tuple,
    corr_wise: str,   
):
    """_summary_."""
    QM = sdba.EmpiricalQuantileMapping.train(
        ref,
        tar,
        nquantiles=quantiles,
        group="time.month",
        kind=corr_wise,
    )
    ds_map = QM.ds
    return ds_map


def historic_correction(
    cl_name: str,
    period: str,
    catalog_fn: Path | str,
    era5_dir: Path | str,
    stime: str,
    etime: str,
    quantiles: list | tuple,
    quantiles_plot: list | tuple,
    output_dir: Path | str,
):
    """_summary_."""
    # Check if directory is there
    output_dir = check_directory(output_dir)
    data_output_dir = check_directory(Path(output_dir, "data", cl_name))
    fig_output_dir = check_directory(Path(output_dir, "figures", cl_name))
    # Set to Path
    catalog_fn = Path(catalog_fn)

    # Get and set the input data
    dc = DataCatalog(data_libs=[catalog_fn])
    ref = xr.open_zarr(Path(era5_dir, f"{cl_name}.zarr"))
    ref = ref.sel(time=slice(stime, etime))
    target = dc.get_rasterdataset(f"{cl_name}_{period}")
    if target.time.dt.calendar == "360_day":
        dtime = datetime.datetime.strptime(etime, "%Y-%m-%dT%H:%M:%S")
        if dtime.day == 31:
            dtime = dtime - datetime.timedelta(days=1)
        etime = dtime.strftime(format="%Y-%m-%dT%H:%M:%S")
    target = target.sel(time=slice(stime, etime))

    # load the catalog yaml
    with open(catalog_fn, "r") as _r:
        catalog = yaml.safe_load(_r)

    cat_entry = copy.deepcopy(CMIP6_CATALOG_BP)
    cat_entry["driver"] = "zarr"
    cat_entry["kwargs"]["chunks"] = "auto" 
    cat_entry["meta"]["info"] = "Bias corrected based on ERA5"

    # Empty outgoing dataset
    out_ds = None
    for var in target:
        # Get the DataArray's
        da_tar = target[var]
        da_ref = ref[var]

        # set correct unit for xclim
        da_tar = da_tar.assign_attrs(
            {"units": UNITS_MAP[var]},
        )
        da_ref = da_ref.assign_attrs(
            {"units": UNITS_MAP[var]},
        )
        da_tar = da_tar.chunk(chunks={"time": len(da_tar.time), "lon": 2, "lat": 2})
        da_ref = da_ref.chunk(chunks={"time": len(da_ref.time), "lon": 2, "lat": 2})

        # Convert the calendar of the reference data
        da_ref = da_ref.convert_calendar(
            da_tar.time.dt.calendar,
            align_on="year",
        )

        # In case of precip, filter out negative values from reference data
        if var in ["precip", "kin"]:
            da_ref = da_ref.where(da_ref > 0, 0)

        # Train the mapping
        corr_wise = UNITS_MULT[var]
        ds_map = train_mapping(
            da_ref,
            da_tar,
            quantiles,
            corr_wise,
        )

        # Backfill all values where nan, inf's are encountered
        ds_map["af"] = ds_map.af.where(np.isfinite(ds_map.af))
        ds_map["af"] = ds_map.af.bfill(dim="quantiles")

        # Bias correct the data
        da_corr = var_correction(da_tar, ds_map=ds_map)
        da_corr = da_corr.transpose("time", "lat", "lon")
        
        # Save the bias correction factors
        ds_map.to_netcdf(
            Path(data_output_dir, f"{var}.nc"),
        )

        # Create figures of the bias correction
        for _q in quantiles_plot:
            fig = plot_spatial_correction(
                da_ref,
                da_tar,
                da_corr,
                var,
                cl_name,
                quantile=_q
            )
            # Set the figure output directory
            out_q_str = str(_q).replace(".", "")
            fig.savefig(
                Path(
                    fig_output_dir,
                    f"{var}_{out_q_str}.png",
                ),
                dpi=300,
            )
        
        # Handle output data
        da_corr.name = var
        if out_ds is None:
            out_ds = da_corr.to_dataset()
            continue
        out_ds[var] = da_corr
    out_ds = out_ds.assign_attrs(target.attrs)

    # Export to archive
    out_fn = Path(output_dir, f"{cl_name}.zarr")
    export_to_zarr(
        out_ds,
        chunks={"time":3000, "lon": len(out_ds.lon), "lat": len(out_ds.lat)},
        out_fn=out_fn,
    )

    cat_path = make_path_relative(out_fn, Path(catalog_fn.parent, "CLIMATE"))
    cat_entry["path"] = cat_path.as_posix()
    catalog[f"{cl_name}_{period}_bc"] = cat_entry

    with open(catalog_fn, "w") as _w:
        yaml.dump(catalog, _w, Dumper=MyDumper, sort_keys=False)
    pass


def future_correction(
    cl_name: str,
    period: str,
    catalog_fn: Path | str,
    ref_dir: Path | str,
    training_data: Path | str,
    stime: str,
    etime: str,
    ref_stime: str,
    ref_etime: str,
    quantiles_plot: tuple | list,
    output_dir: Path | str,
):
    """_summary_."""
    # Check if directory is there
    output_dir = check_directory(output_dir)
    fig_output_dir = check_directory(Path(output_dir, "figures", cl_name))
    # Set to Path
    catalog_fn = Path(catalog_fn)

    # Get and set the input data
    dc = DataCatalog(data_libs=[catalog_fn])
    ref = xr.open_zarr(Path(ref_dir, f"{cl_name}.zarr"))
    ref = ref.sel(time=slice(ref_stime, ref_etime))
    target = dc.get_rasterdataset(f"{cl_name}_{period}")
    if target.time.dt.calendar == "360_day":
        dtime = datetime.datetime.strptime(etime, "%Y-%m-%dT%H:%M:%S")
        if dtime.day == 31:
            dtime = dtime - datetime.timedelta(days=1)
        etime = dtime.strftime(format="%Y-%m-%dT%H:%M:%S")
    target = target.sel(time=slice(stime, etime))
            
    # load the catalog yaml
    with open(catalog_fn, "r") as _r:
        catalog = yaml.safe_load(_r)

    cat_entry = copy.deepcopy(CMIP6_CATALOG_BP)
    cat_entry["driver"] = "zarr"
    cat_entry["kwargs"]["chunks"] = "auto"
    cat_entry["meta"]["info"] = "Bias corrected based on ERA5"

    # Empty outgoing dataset
    out_ds = None
    for var in target:
        # Get the DataArray's
        da_tar = target[var]
        da_ref = ref[var]

        # set correct unit for xclim
        da_tar = da_tar.assign_attrs(
            {"units": UNITS_MAP[var]},
        )
        da_ref = da_ref.assign_attrs(
            {"units": UNITS_MAP[var]},
        )
        da_tar = da_tar.chunk(chunks={"time": len(da_tar.time), "lon": 2, "lat": 2})
        da_ref = da_ref.chunk(chunks={"time": len(da_ref.time), "lon": 2, "lat": 2})

        # Get the trained data
        ds_map = xr.open_dataset(
            Path(training_data, cl_name, f"{var}.nc")
        )
        ds_map = ds_map.chunk(chunks={"lon": 2, "lat": 2})

        # Bias correct the data
        da_corr = var_correction(da_tar, ds_map=ds_map)
        da_corr = da_corr.transpose("time", "lat", "lon")
        
        # Create figures of the bias correction
        for _q in quantiles_plot:
            fig = plot_spatial_correction(
                da_ref,
                da_tar,
                da_corr,
                var,
                cl_name,
                quantile=_q
            )
            # Set the figure output directory
            out_q_str = str(_q).replace(".", "")
            fig.savefig(
                Path(
                    fig_output_dir,
                    f"{var}_{out_q_str}.png",
                ),
                dpi=300,
            )
        
        # Handle output data
        da_corr.name = var
        if out_ds is None:
            out_ds = da_corr.to_dataset()
            continue
        out_ds[var] = da_corr
    out_ds = out_ds.assign_attrs(target.attrs)

    # Export to archive
    out_fn = Path(output_dir, f"{cl_name}.zarr")
    export_to_zarr(
        out_ds,
        chunks={"time":3000, "lon": len(out_ds.lon), "lat": len(out_ds.lat)},
        out_fn=out_fn,
    )

    cat_path = make_path_relative(out_fn, Path(catalog_fn.parent, "CLIMATE"))
    cat_entry["path"] = cat_path.as_posix()
    catalog[f"{cl_name}_{period}_bc"] = cat_entry

    with open(catalog_fn, "w") as _w:
        yaml.dump(catalog, _w, Dumper=MyDumper, sort_keys=False)


if __name__ == "__main__":
    stime = "1970-01-01T03:00:00"
    etime = "2014-12-31T21:00:00"
    quantiles = np.linspace(0.01,0.99,30)

    for cl_name in CMIP6_MODELS:
        period = "historic"
        hist_output_dir = Path(
            "/p/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE",
            f"{period.upper()}_BC"
        )
        # Historic period
        historic_correction(
            cl_name,
            period=period,
            catalog_fn="/p/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml",
            era5_dir="/p/1000365-002-wflow/tmp/usgs_wflow/data/TMP/era5_climate",
            stime=stime,
            etime=etime,
            quantiles=quantiles,
            quantiles_plot=[0.5, 0.9],
            output_dir=hist_output_dir,
        )
        # Future period
        period="future"
        fut_output_dir = Path(
            "/p/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE",
            f"{period.upper()}_BC"
        )
        stime_fut = "2015-01-01T03:00:00"
        etime_fut = "2050-12-31T21:00:00"
        future_correction(
            cl_name,
            period=period,
            catalog_fn="/p/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml",
            ref_dir="/p/1000365-002-wflow/tmp/usgs_wflow/data/TMP/era5_climate",
            training_data=Path(hist_output_dir, "data"), 
            stime=stime_fut,
            etime=etime_fut,
            ref_stime=stime,
            ref_etime=etime,
            quantiles_plot=[0.5, 0.9],
            output_dir=fut_output_dir,
        )       
    pass