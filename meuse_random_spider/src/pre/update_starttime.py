from hydromt_wflow import WflowModel
import xarray as xr
import pandas as pd


def update_starttime(starttime, dt, calendar):
    # Create timeseries
    dates = xr.cftime_range(
        start=starttime,
        periods=1,
        calendar=calendar,
    )
    # Determine offset
    offset = pd.Timedelta(dt, unit="s")

    # Subtract from first timestep
    return dates[0] - offset


def update_config(root, cmip_model, config_fn):
    # Read model
    mod = WflowModel(
        root=root,
        config_fn=config_fn,
        mode="r+",
    )

    # Update starttime by subtracting one dt
    newtime = update_starttime(
        starttime=mod.config["starttime"],
        dt=mod.config["timestepsecs"],
        calendar=mod.config["calendar"],
    )

    # Update config and overwrite existing file
    mod.config["starttime"] = newtime.strftime("%Y-%m-%dT%H:%M:%S")
    mod.write_config()

    return mod


if __name__ == "__main__":

    basins = ["KING", "PIERCE"]
    cmip_models = ["CNRM", "EcEarth", "GFDL", "HadGemHH", "HadGemHM", "HadGemHMsst"]

    model_dir = R"p:\1000365-002-wflow\tmp\usgs_wflow\models\MODELDATA_{basin}_CLIMATE"

    for basin in basins:
        for cmip_model in cmip_models:

            root = model_dir.format(basin=basin)

            # mod = update_config(
            #     root=root,
            #     cmip_model=cmip_model,
            #     config_fn=f"wflow_sbm_cmip6_{cmip_model}_historic.toml",
            # )

            # mod = update_config(
            #     root=root,
            #     cmip_model=cmip_model,
            #     config_fn=f"wflow_sbm_cmip6_{cmip_model}_future.toml",
            # )

            # mod = update_config(
            #     root=root,
            #     cmip_model=cmip_model,
            #     config_fn=f"wflow_sbm_cmip6_{cmip_model}_historic_bc.toml",
            # )

            # mod = update_config(
            #     root=root,
            #     cmip_model=cmip_model,
            #     config_fn=f"wflow_sbm_cmip6_{cmip_model}_future_bc.toml",
            # )
