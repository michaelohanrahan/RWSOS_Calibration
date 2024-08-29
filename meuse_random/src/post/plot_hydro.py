import hydromt
from hydromt_wflow import WflowModel
import xarray as xr
import numpy as np
import os
import xarray as xr
import pandas as pd
from datetime import datetime
from tqdm import tqdm

from func_plot_signature import plot_signatures, plot_hydro
from func_io import read_filename_txt, read_lakefile


for basin in ["PIERCE", "KING"]:
    main_dir = R"p:\1000365-002-wflow\tmp\usgs_wflow"

    fn_obs = R"p:\1000365-002-wflow\tmp\usgs_wflow\data\GAUGES\discharge_obs.nc"
    output_dir = Rf"{main_dir}\models\MODELDATA_{basin}_500M\run_default"
    figure_dir = Rf"{main_dir}\docs\images\calib\{basin}"
    ksat_values = [
        5,
        25,
        100,
        500,
        1000,
    ]
    # 5, 10, 20, 25, 50, 75, 100, 200, 250, 500, 750, 1000]

    start_long, end_long = "2012-01-01", "2018-12-31"
    start_1, end_1 = "2012-09-01", "2013-08-31"
    start_2, end_2 = "2013-09-01", "2014-08-31"
    start_3, end_3 = "2015-09-01", "2016-08-31"

    model_runs = dict()
    for ksat_value in ksat_values:
        model_runs[f"ksat_{ksat_value}"] = f"{output_dir}/ksat_{ksat_value}/output.csv"

    idx_start = 365 * 2 * int(24 / 3)

    toml_default_fn = f"{output_dir}/../blueprint_wflow_sbm.toml"

    gauges_maps = [
        "gauges_obs",
    ]

    # Get stations within model
    root = os.path.dirname(toml_default_fn)
    mod = WflowModel(root, config_fn=os.path.basename(toml_default_fn), mode="r")
    # mod_stations = mod.staticgeoms["gauges_gauges-obs_hr"].stations.values

    mod_stations = np.array([], dtype=int)
    for gauge_map in gauges_maps:
        mod_stations = np.append(
            mod_stations, mod.geoms[gauge_map]["SiteNumber"].values
        )

    # # Get names of the locations
    # fn_names = R"p:\11205237-grade\wflow\wflow_rhine_julia\measurements\vanBart\discharge_obs_hr_appended_station_list.csv"
    # df_locs = pd.read_csv(fn_names, index_col=0, encoding="unicode_escape")
    # df_locs = df_locs.loc[mod_stations, :]
    # df_locs["names"] = df_locs.station_names.str.split("'", expand=True)[1]
    # # Convert to dictionary
    # stations_dict = df_locs.names.to_dict()

    # Read observations
    ds_obs = xr.open_dataset(fn_obs)
    df_obs_tmp = ds_obs.to_dataframe()
    time_index = np.unique(df_obs_tmp.index.get_level_values("time"))

    # Get names of the stations
    stations_dict = dict()
    for station in ds_obs.data_vars:
        stations_dict[station] = dict()
        name = ds_obs[station].attrs["site_name"].title()
        name = name.replace(", Wa", ", WA")
        stations_dict[station]["name"] = name

        area = ds_obs[station].attrs["drainage_area"]
        stations_dict[station]["area"] = area

    # Convert to DataFrame
    df_obs = pd.DataFrame(index=df_obs_tmp.index)
    for station in mod_stations:
        try:
            df_obs[f"Q_{station}"] = df_obs_tmp[str(station)] / 35.315
        except:
            df_obs[f"Q_{station}"] = np.nan

    ### prepare dataset to make plots
    colors = [
        "#a6cee3",
        "#1f78b4",
        "#b2df8a",
        "#33a02c",
        "orange",
        "#fb9a99",
        "#e31a1c",
        "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#ffff99",
        "#b15928",
    ]

    # save output in a dictionary
    runs_dict = {}
    for key in model_runs.keys():
        runs_dict[key] = pd.read_csv(
            model_runs[key], index_col=0, header=0, parse_dates=True
        )[idx_start:]

    ### #make dataset
    variables = ["Q"]
    runs = ["Obs.", *model_runs.keys()]
    rng = runs_dict[key].index
    rng.freq = "3H"

    S = np.zeros((len(rng), len(mod_stations), len(runs)))
    v = (("time", "stations", "runs"), S)
    h = {k: v for k in variables}

    ds = xr.Dataset(
        data_vars=h,
        coords={"time": rng, "stations": mod_stations, "runs": runs},
    )
    ds = ds * np.nan

    # fill dataset with model and observed data
    ds["Q"].loc[dict(runs="Obs.")] = df_obs.loc[rng]
    for key, item in runs_dict.items():
        for sub in list(map(str, mod_stations)):
            ds["Q"].loc[dict(runs=key, stations=int(sub))] = item.loc[rng, f"Q_{sub}"]

    ##%%
    plot_colors = colors[: len(runs_dict)]
    for station_id in tqdm(mod_stations):
        # Extract station name
        try:
            tmp_name = stations_dict[str(station_id)]["name"]
            tmp_area = stations_dict[str(station_id)]["area"]
            station_name = f"{tmp_name} ({tmp_area} mi$^2$)"
        except KeyError:
            station_name = str(station_id)

        # Get data for that stations
        dsq = ds.sel(stations=station_id)
        # Prep labels
        labels = np.delete(dsq.runs.values, np.where(dsq.runs.values == "Obs."))

        # plot hydro

        plot_hydro(
            dsq=dsq,
            start_long=start_long,
            end_long=end_long,
            start_1=start_1,
            end_1=end_1,
            start_2=start_2,
            end_2=end_2,
            start_3=start_3,
            end_3=end_3,
            labels=labels,
            colors=plot_colors,
            Folder_out=figure_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
        )

        plot_signatures(
            dsq=dsq,
            labels=labels,
            colors=plot_colors,
            Folder_out=figure_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
        )
