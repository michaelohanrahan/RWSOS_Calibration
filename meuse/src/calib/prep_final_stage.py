import shutil
from pathlib import Path

import xarray as xr

from set_config import main as create_config


def main(
    in_perf: tuple | list,
    in_cfg: Path | str,
    in_maps: Path | str,
    cfg_args: tuple | list,
    out_cfg: Path | str,
    out_maps: Path | str,
    out_perf: Path | str,
):
    """_summary_"""
    # Copy the staticmaps
    out_dir = Path(out_maps).parent

    if not out_dir.exists():
        out_dir.mkdir()

    shutil.copy2(
        in_maps,
        out_maps,
    )

    # Set the config file
    cfg_args += [[out_cfg],]

    create_config(
        in_cfg,
        *cfg_args,
    )

    out_ds = None

    for _ds in in_perf:
        temp_ds = xr.open_dataset(_ds)
        if out_ds is None:
            out_ds = temp_ds.copy()
            continue
        out_ds = xr.concat([out_ds, temp_ds], dim="gauges")

    out_ds.to_netcdf(out_perf)


if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        
        main(
            mod.input.performance,
            mod.params.cfg_template,
            mod.params.staticmaps,
            mod.params.cfg_args,
            mod.output.cfg,
            mod.output.staticmaps,
            mod.output.performance,
        )
    
    else:
        pass