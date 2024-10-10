import calendar
import copy
import datetime
import re
import yaml
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from fiat.log import setup_default_log
from hydromt import DataCatalog
from hydromt.gis_utils import meridian_offset

from puget import (
    CMIP6_CATALOG_BP,
    CMIP6_UNIT_ADD,
    CMIP6_UNIT_MULT, 
    CMIP6_VARMAP,
    MyDumper,
    check_directory, 
    export_to_zarr,
)


def create_yearly_cycle(
    grouped: pd.DataFrame,
    year: int,
):
    """_summary_"""
    cycle = []
    for m in range(1, 13):
        _, nod = calendar.monthrange(year, m)
        m_vals = grouped[grouped.index.month==m]
        m_not = f"{m:02d}"
        df = pd.concat(
            [m_vals,]*nod
        )
        df.index = pd.date_range(
            f"{year}-{m_not}-01T00:00:00",
            f"{year}-{m_not}-{nod}T21:00:00",
            freq="3H",
        )
        cycle.append(df)
    
    cycle = pd.concat(cycle)
    return cycle


def zonal_stats(
    ds_name: str,
    region: tuple | list,
):
    """_summary_."""

    d = DataCatalog(
        data_libs=["/p/1000365-002-wflow/tmp/usgs_wflow/data/deltares_data.yml"]
    )
    era5 = d.get_rasterdataset(
        ds_name,
        time_tuple=("1990-01-01T01:00:00", "2020-01-01T21:00:00")
    )
    
    kin = meridian_offset(era5.kin)
    # Shift according to timezone
    # kin["time"] = kin["time"] - pd.Timedelta(hours=8)
    # Get the region for zonal stats
    region = [gpd.read_file(item).dissolve() for item in region]
    region = pd.concat(region)
    region = region.dissolve()
    
    # Create a 1D timeseries from the region
    # Sample this already to a 3hourly temporal resolution
    zonal = kin.raster.zonal_stats(region, "mean")
    zonal = zonal.sel(index=0).drop("index")
    zonal = zonal.to_dataframe()
    zonal = zonal.resample("3H", label="right", closed="right").sum()
    zonal.drop("spatial_ref", inplace=True, axis=1)
    zonal[zonal<0] = 0

    # Group the data
    grouped = zonal.groupby([zonal.index.month, zonal.index.hour]).mean()
    grouped.index.names = ["month", "hour"]
    monthly_sum = grouped.groupby(level=[0]).sum()
    grouped = grouped / monthly_sum
    # Resize it to a yearly series
    grouped = grouped.reset_index()
    grouped['time'] = grouped.apply(lambda row: pd.Timestamp(year=2000, month=int(row['month']), day=1, hour=int(row['hour'])), axis=1)
    grouped.set_index("time", inplace=True)
    grouped.drop(["month", "hour"], axis=1, inplace=True)
    grouped.rename({"kin_mean": "kin"}, inplace=True, axis=1)

    # non_leap = create_yearly_cycle(grouped, 2001)
    # leap = create_yearly_cycle(grouped, 2000)

    patterns = {
        _m: grouped[grouped.index.month==_m].values[:,0].tolist() 
        for _m in range(1,13)
    }
    for key, item in patterns.items():
        m_pattern = sum(item)/len(item)
        item = [_i/m_pattern for _i in item]
        patterns[key] = item
    return patterns #, non_leap, leap


def resample_daily(
    ds: xr.Dataset,
    patterns: pd.DataFrame,
    var: str,
):
    """_summary_."""

    total = []
    lon_len = ds.lon.__len__()
    lat_len = ds.lat.__len__()

    for _m in range(1,13):
        pattern = patterns[_m]
        # Upscale the data
        monthly_vals = ds.sel(time=ds.time.dt.month == _m)
        len_monthly = monthly_vals.time.__len__()
        monthly_vals = monthly_vals.resample(time="3H").asfreq()
        monthly_vals = monthly_vals[var]
        # Add the dates in front and back
        first_steps = [
            monthly_vals['time'].values[0] - datetime.timedelta(hours=hour)
            for hour in [12,9,6,3]
        ]
        da_first = xr.DataArray(dims=["time"], coords={"time": first_steps})
        back_steps = [
            monthly_vals['time'].values[-1] + datetime.timedelta(hours=hour)
            for hour in [3,6,9]
        ]
        da_last = xr.DataArray(dims=["time"], coords={"time": back_steps})
        monthly_vals = xr.concat([da_first, monthly_vals, da_last], dim="time")
        # Back and forward fill (4 timesteps back, 3 forward)
        monthly_vals = monthly_vals.bfill(dim="time", limit=4)
        monthly_vals = monthly_vals.ffill(dim="time", limit=3)
        monthly_vals = monthly_vals.dropna(dim="time")
        
        d3 = np.array(
            [[[item for _ in range(lon_len)] for _ in range(lat_len)] 
             for item in pattern]
        )
        mult = np.tile(d3, (len_monthly, 1, 1))
        monthly_vals *= mult
        monthly_vals.name = var
        total.append(monthly_vals.to_dataset())
        del monthly_vals
        del mult

    total = xr.concat(total, dim="time")
    total[var] = total[var].assign_attrs(ds[var].attrs)
    return total


def reset_calendar(
    ds: xr.Dataset,
    starttime: str,
    deltat: int,
    length: int,
):
    """_summary_."""
    dates = xr.cftime_range(
        start=starttime,
        periods=length,
        freq=f"{deltat}H",
        calendar="360_day",
    )

    # Assign it
    ds = ds.assign_coords(
        {"time": dates},
    )
    return ds


def process_data(
    dir_fn: Path| str,
    period: str,
    pattern: str,
    out_dir: Path | str,
    out_catalog: Path | str,
):
    """_summary_."""
    logger.info(f"Setting up {period} climate data processing")

    # Setup the directories and pattern
    dir_fn = Path(dir_fn)
    out_dir = Path(out_dir)
    pattern = re.compile(pattern)

    # Setup the catalog
    logger.info("Create or load the climate data catalog.")
    catalog = {}
    catalog_meta = {
        "meta": {
            "root": out_dir.parent.as_posix(),
            "version": 2024.1,
        },
    }
    out_catalog = Path(out_catalog)
    if out_catalog.exists():
        with open(out_catalog, "r") as _r:
            catalog = yaml.safe_load(_r)
    if "meta" not in catalog:
        catalog.update(catalog_meta)

    logger.info("Calculating monthly patterns for daily resampling")
    # Monthly patterns
    patterns = zonal_stats(
        "era5_hourly_zarr",
        [
            "/p/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_500M/staticgeoms/basins.geojson",
            "/p/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_PIERCE_500M/staticgeoms/basins.geojson"
        ],
    )

    for _f in dir_fn.iterdir():
        matches = pattern.match(_f.name)
        if matches is None:
            continue

        ds_name, var = matches.groups()
        if var not in CMIP6_VARMAP:
            continue

        logger.info(f"Processing model: {ds_name}; variable: {var}")

        ds = xr.open_dataset(_f)
        out_subdir = check_directory(Path(out_dir, ds_name))

        output_name = f"{ds_name}_{period}"

        if ds.lon.values[0] > 180:
            ds = ds.assign_coords({"lon": ds.lon.values - 360})
        if ds.lat.values[-1] > ds.lat.values[0]:
            ds = ds.raster.flipud()

        dt = (ds.time.values[1] - ds.time.values[0])
        dt = np.timedelta64(dt, "h").astype(int)

        logger.info(f"Timestep: {dt} hours")

        # Shift the time axis if timestep fall on 30 min
        if ds.time.dtype.char == "M":
            sd = pd.Timestamp(ds.time.values[0]).strftime("%Y-%m-%dT%H:%M:%S")
            dates = xr.cftime_range(
                sd, 
                periods=len(ds.time), 
                freq=f"{dt}H", 
                calendar="proleptic_gregorian",
            )
            ds = ds.assign_coords(
                {"time": dates},
            )

        sd = ds.time.values[0]
        ed = ds.time.values[-1]

        logger.info(f"Starttime: {sd}")

        nyear = ed.year - sd.year
        if ed.month > 10:
            nyear += 1
        nod = len(ds.time) / (24/dt) / nyear

        logger.info(f"Number of days per year: {nod}")

        if sd.minute != 0:
            ds['time'] = ds['time'] + datetime.timedelta(hours=1.5)
            logger.warning(f"Readjusted starttime: {ds.time.values[0]}")

        # Reset time axis to proper 360 day calendar if there are 360 days in one year
        if nod < 361:
            logger.warning("Calendar detected with 360 days, \
settings months with 30 days")
            sd = ds.time.values[0]
            ed = ds.time.values[-1]
            logger.info(f"Original calendar range {sd} - {ed}")
            ds = reset_calendar(
                ds,
                sd.strftime("%Y-%m-%dT%H:%M:%S"),
                dt,
                len(ds.time),
            )
            sd = ds.time.values[0]
            ed = ds.time.values[-1]
            logger.info(f"New 360 calendar range {sd} - {ed}")

        # Resample if not at a 3 hourly resolution
        if dt == 6:
            logger.warning("Resamling from 6 hourly to 3 hourly by backfilling")
            first_time_step = ds['time'].values[0] - datetime.timedelta(hours=3)
            ds_first = xr.Dataset({}, coords={'time': [first_time_step]})

            ds = xr.concat([ds_first, ds], dim='time')
            ds[var] = ds[var].bfill(dim="time", limit=1)
            ds = ds.resample(time='3H').bfill()
        elif dt == 24:
            logger.warning("Resamling from daily to 3 hourly based on ERA5")
            ds = resample_daily(
                ds,
                patterns,
                var,
            )

        if output_name not in catalog: 
            catalog[output_name] = copy.deepcopy(CMIP6_CATALOG_BP)
        catalog[output_name]["path"] = Path(
            period.upper(), 
            ds_name, 
            "*.nc"
        ).as_posix()
        
        hydromt_var = CMIP6_VARMAP[var]
        if "rename" not in catalog[output_name]:
            catalog[output_name]["rename"] = {}
        catalog[output_name]["rename"][var] = hydromt_var
        if "unit_add" not in catalog[output_name]:
            catalog[output_name]["unit_add"] = {}
        if hydromt_var in CMIP6_UNIT_ADD:
            catalog[output_name]["unit_add"][hydromt_var] = CMIP6_UNIT_ADD[
                hydromt_var
            ]
        if "unit_mult" not in catalog[output_name]:
            catalog[output_name]["unit_mult"] = {}
        if hydromt_var in CMIP6_UNIT_MULT:
            catalog[output_name]["unit_mult"][hydromt_var] = CMIP6_UNIT_MULT[
                hydromt_var
            ]

        out_file_path = Path(out_subdir, f"{var}.nc")
        logger.info(f"Writing file to {out_file_path.as_posix()}")
        ds.to_netcdf(
            out_file_path,
            encoding={var: {"zlib": True, "dtype": "float64"}},
        )
    
    logger.info(
        f"Writing blueprint catalog to f{out_catalog.as_posix()}"
    )

    with open(out_catalog, "w") as _w:
        yaml.dump(catalog, _w, Dumper=MyDumper, sort_keys=False)

    logger.info(f"Done with processing of {period} climate data!")


def zarr_achives(
    p: Path | str,
    c: Path | str,
    period: str,
):
    """_summary_"""
    logger.info("Converting processed netCDF data to zarr archives")
    p = Path(p)

    c_path = Path(c)
    with open(c_path, "r") as _r:
        catalog = yaml.safe_load(_r)

    for _p in p.iterdir():
        p_name = _p.stem
        if not _p.is_dir():
            continue

        logger.info(f"Converting {p_name}")

        catalog_name = f"{p_name}_{period}"

        ds = xr.open_mfdataset(
            Path(_p, "*.nc").as_posix(),
        )

        out_file_path = Path(p, f"{p_name}.zarr")
        logger.info(f"Writing archive to {out_file_path.as_posix()}")
        export_to_zarr(
            ds,
            chunks={"time": 3000, "lon": len(ds.lon), "lat": len(ds.lat)},
            out_fn=out_file_path.as_posix(),
        )

        catalog[catalog_name]["path"] = Path(
            period.upper(),
            f"{p_name}.zarr",
        ).as_posix()
        catalog[catalog_name]["driver"] = "zarr"
        catalog[catalog_name]["kwargs"]["chunks"] = "auto"

    logger.info("Adjusting the blueprint catalog")
    with open(c_path, "w") as _w:
        yaml.dump(catalog, _w, Dumper=MyDumper, sort_keys=False)

    logger.info("Done converting to zarr!")


if __name__ == "__main__":
    period = "future"
    bp_catalog = "/p/1000365-002-wflow/tmp/usgs_wflow/data/climate_catalog.yml"

    output_dir = check_directory(
        f"/p/1000365-002-wflow/tmp/usgs_wflow/data/CLIMATE/{period.upper()}"
    )

    logger = setup_default_log(
        "climate_processing",
        2,
        output_dir
    )

    process_data(
        f"/p/11208413-usgs-coop-22-23/PugetSound/data/CMIP6/hybrid/{period}",
        period,
        f"(.*)_{period}_(.*).nc",
        output_dir,
        bp_catalog,
    )

    zarr_achives(
        output_dir,
        bp_catalog,
        period,
    )