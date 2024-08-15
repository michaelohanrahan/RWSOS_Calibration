import xarray as xr
from pathlib import Path
from setuplog import setup_logging
import traceback

def main(
    l,
    files: list,
    out: str|Path,
):
    l.info(f"Combining {len(files)} files to {out}")
    ds = xr.open_mfdataset(files, combine="by_coords")
    ds.to_netcdf(out)
    with open(out.parent / "eval.done", "w") as f:
        f.write("")
    l.info(f"Saved {out}")


if __name__ == "__main__":
    l = setup_logging("data/0-log", "05-combine_evaluated.log")
    try:
        if 'snakemake' in globals():
            main(
                l,
                snakemake.input.files,
                snakemake.output,
            )
        else:
            main(l,
                [
                    r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\ksat~0.5\f~0.75\rd~1.0\st~0.5\nr~0.5\ml~0.3\evaluated.nc",
                    r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\ksat~0.75\f~0.5\rd~1.75\st~1.0\nr~1.0\ml~0.3\evaluated.nc",
                ],
                out=r"p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\data\2-interim\calib_data\level0\performance.nc",
            )
    except Exception as e:
        l.exception(e)
        l.error(traceback.format_exc())
        raise