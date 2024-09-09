from pathlib import Path

import numpy as np
import xarray as xr

from func_plot_signature import plot_signatures, plot_hydro


def main(
    md_data: Path | str,
    obs_data: Path | str,
    gauges: Path | str,
    starttime: Path | str,
    endtime: Path | str,
    period_start: tuple | list,
    period_length: tuple | list,
    period_unit: str,
    output_dir: Path | str,
):
    """_summary_"""
    # cfg path
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir()

    md_ds = xr.open_dataset(md_data)
    obs_ds = xr.open_dataset(obs_data)
    obs_ds.Q.values = obs_ds.Q.values * (0.3048**3)

    temp = obs_ds.sel(Q_gauges_obs=gauges, time=slice(starttime, endtime))
    da = xr.DataArray(
        np.expand_dims(temp.Q.values, 2),
        coords = {
            "time": temp.time.values,
            "stations": gauges,
            "runs": ["Obs.",]
        },
        dims = ["time", "stations", "runs"]
    )
    temp = md_ds.sel(Q_gauges_obs=gauges, time=slice(starttime, endtime))
    da = xr.concat(
        [
            da, 
            xr.DataArray(
                np.expand_dims(temp.Q.values, 2),
                coords = {
                    "time": temp.time.values,
                    "stations": gauges,
                    "runs": ["Md.",]
                },
                dims = ["time", "stations", "runs"],
            ),
        ],
        dim="runs",
    )
    ds = da.to_dataset(name="Q")

    hydro_periods = []
    for idx, _s in enumerate(period_start):
        _s = np.datetime64(_s)
        _e = _s + np.timedelta64(period_length[idx], period_unit)
        hydro_periods.append(
            (
                np.datetime_as_string(_s, unit="D"),
                np.datetime_as_string(_e, unit="D"),
            )
        )
    
    pass
    ### prepare dataset to make plots
    colors = [
        # "#a6cee3",
        # "#1f78b4",
        # "#b2df8a",
        # "#33a02c",
        # "orange",
        # "#fb9a99",
        # "#e31a1c",
        # "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#ffff99",
        "#b15928",
    ]

    # ##%%
    plot_colors = colors[:1]
    for station_id in gauges:
        # Extract station name
        try:
            station_name = str(station_id)
        except KeyError:
            station_name = str(station_id)

        # Get data for that stations
        dsq = ds.sel(stations=station_id)
        # Prep labels
        labels = np.delete(dsq.runs.values, np.where(dsq.runs.values == "Obs."))

        # Plot hydrograpgh
        plot_hydro(
            dsq=dsq,
            start_long=starttime.split("T")[0],
            end_long=endtime.split("T")[0],
            start_1=hydro_periods[0][0],
            end_1=hydro_periods[0][1],
            start_2=hydro_periods[1][0],
            end_2=hydro_periods[1][1],
            start_3=hydro_periods[2][0],
            end_3=hydro_periods[2][1],
            labels=labels,
            colors=plot_colors,
            Folder_out=output_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
        )

        # plot signatures
        plot_signatures(
            dsq=dsq,
            labels=labels,
            colors=plot_colors,
            Folder_out=output_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
        )


if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        main(
            mod.input.scalar,
            mod.params.observed_data,
            mod.params.gauges,
            mod.params.starttime,
            mod.params.endtime,
            mod.params.period_startdate,
            mod.params.period_length,
            mod.params.period_unit,
            mod.params.output_dir,
        )
    
    else:
        main(
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/Intermittent_result/run_default/output_scalar.nc",
            "p:/1000365-002-wflow/tmp/usgs_wflow/data/GAUGES/discharge_obs_combined.nc",
            ['12112600', '12108500', '12105900', '12117000', '12115500', '12114500', '12120600', '12120000', '12106700', '12115000', '12121600'],
            "2012-01-01T00:00:00",
            "2018-12-31T00:00:00",
            ["2012-09-01", "2013-09-01", "2015-09-01"],
            [365, 365, 365],
            "D",
            "p:/1000365-002-wflow/tmp/usgs_wflow/models/MODELDATA_KING_CALIB/Intermittent_result/figures"
        )
