import copy
from pathlib import Path

import tomli
import tomli_w
from setuplog import setup_logging
import traceback


def main(
    level: int,
    cfg: Path | str,
    starttime: str,
    endtime: str,
    # timestep: str | int,
    forcing_path: Path | str,
    out: tuple | list,
    gaugemap: str,
    l
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
    l.info(f"Loading blueprint configuration from {cfg}")
    with open(cfg, "rb") as _r:
        data = tomli.load(_r)

    for out_file in out:
        l.info(f"target dir: {out_file}")
        
        out_cfg = copy.deepcopy(data)

        # Ensure the directory
        out_file_dir = Path(out_file).parent
        if not out_file_dir.exists():
            out_file_dir.mkdir()
        
        out_cfg["starttime"] = starttime
        out_cfg["endtime"] = endtime
        # out_cfg["timestepsecs"] = timestep

        # Set the forcing path to the source dir
        out_cfg["input"]["path_forcing"] = forcing_path.as_posix()
        
        #ex: p:/11209265-grade2023/wflow/RWSOS_Calibration/meuse/data/2-interim/calib_data/level0/ksat~0.4/f~1.5/st~1.5/nr~0.5/rd~0.8/ml~0.0/nl~0.8/wflow_sbm.toml
        # level = int(out_file.split('/')[-8][-1]) # the integer of the level
        thickness = str(''.join(out_file.split('/')[-7].split('~')[-1].split('.'))) #the joined float value as a string
        l.info(f"level: {level}, thickness: {thickness}")
        
        #ex: "instate_level{level}_ST{thickness}.nc"
        st_instate_path = Path(forcing_path).parent / "instates" / f"instate_level{level}_ST{thickness}.nc"
        
        #change the instate depending on the soilthickness
        #if instate is not an existing key, add it
        if "path_input" not in out_cfg["state"]:
            out_cfg["state"]["path_input"] = {}
            
        out_cfg["state"]["path_input"] = st_instate_path.as_posix()
        out_cfg["input"]["path_static"] = "staticmaps.nc"
        
        out_cfg.pop('csv')
        
        if "netcdf" not in out_cfg:
            out_cfg["netcdf"] = {}
        
        out_cfg["netcdf"]["path"] = "output_scalar.nc"
        out_cfg["netcdf"]["variable"] = [
            {
                "name":"Q",
                "map":gaugemap,
                "parameter":"lateral.river.q_av"
            },
            {
                "name":"Q_hbv",
                "map":"gauges_hbv",
                "parameter":"lateral.river.q_av"
            }
        ]
        
        l.info(f"st_instate_path: {st_instate_path.as_posix()}")
        l.info(f"writing to {out_file}")
        

        # Write the settings file
        with open(out_file, "wb") as _w:
            tomli_w.dump(out_cfg, _w)
        
        assert Path(out_file).exists(), f"Failed to write the configuration file to {out_file}"

if __name__ == "__main__":
    l = setup_logging("data/0-log", "02-set_config.log")
    if "snakemake" in globals():
        mod = globals()["snakemake"]
        try:
            main(
                mod.params.level,
                mod.params.cfg_template,
                mod.params.starttime,
                mod.params.endtime,
                # mod.params.timestep,
                mod.params.forcing_path,
                mod.output,
                mod.params.gaugemap,
                l
            )
        except Exception as e:
            l.error(f"An error occurred: {e}")
            l.error(traceback.format_exc())
            raise e
    else:
        pass

