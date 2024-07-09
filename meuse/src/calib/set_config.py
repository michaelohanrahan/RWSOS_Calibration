import copy
from pathlib import Path

import tomli
import tomli_w


def main(
    cfg: Path | str,
    starttime: str,
    endtime: str,
    timestep: str | int,
    forcing_path: Path | str,
    out: tuple | list,
):
    """
    This script modifies a blueprint configuration file for a specific time period and forcing path.

    It takes the following arguments:
    - `cfg`: The path to the blueprint configuration file.
    - `starttime`: The start time for the simulation.
    - `endtime`: The end time for the simulation.
    - `timestep`: The time step for the simulation.
    - `forcing_path`: The path to the forcing data.
    - `out`: A list of paths where the modified configuration files will be written.

    The script performs the following steps for each path in `out`:
    1. Loads the blueprint configuration file.
    2. Creates a deep copy of the configuration.
    3. Sets the start time, end time, time step, and forcing path in the copied configuration.
    4. Writes the modified configuration to the output file path.

    The script ensures that the directory for each output file path exists before writing the file. If the directory does not exist, it is created.
    """
    # Load the blueprint
    with open(cfg, "rb") as _r:
        data = tomli.load(_r)

    for out_file in out:

        out_cfg = copy.deepcopy(data)

        # Ensure the directory
        out_file_dir = Path(out_file).parent
        if not out_file_dir.exists():
            out_file_dir.mkdir()
        
        out_cfg["starttime"] = starttime
        out_cfg["endtime"] = endtime
        out_cfg["timestepsecs"] = timestep

        # Set the forcing path to the source dir
        out_cfg["input"]["path_forcing"] = forcing_path.as_posix()

        # Write the settings file
        with open(out_file, "wb") as _w:
            tomli_w.dump(out_cfg, _w)

if __name__ == "__main__":
    if "snakemake" in globals():
        mod = globals()["snakemake"]

        main(
            mod.params.cfg_template,
            mod.params.starttime,
            mod.params.endtime,
            mod.params.timestep,
            mod.params.forcing_path,
            mod.output,
        )

    else:
        pass

