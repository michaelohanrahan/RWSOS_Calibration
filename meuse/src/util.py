"""Some utility."""
from os.path import relpath
from pathlib import Path, WindowsPath

import yaml

CMIP6_CATALOG_BP = {
    "crs": 4326,
    "data_type": "RasterDataset",
    "driver": "raster",
    "kwargs": {"chunks": {"time": 3000}},
    "meta": {"category": "Climate scenario", "source_author": "USGS"},
    "path": "",
}

CMIP6_MODELS = [
    "cmip6_CNRM",
    "cmip6_EcEarth",
    "cmip6_GFDL",
    "cmip6_HadGemHH",
    "cmip6_HadGemHM",
    "cmip6_HadGemHMsst",
]

CMIP6_UNIT_ADD = {
    "temp": -273.15,
}

CMIP6_UNIT_MULT = {
    "kin": 0.0001,
    "precip": 10800,
    "press_msl": 1.0e-6,
    "temp": 0.0001,
}

CMIP6_VARMAP = {
    "pr": "precip",
    "psl": "press_msl",
    "rsds": "kin",
    "tas": "temp",
}

UNITS_MAP = {
    "precip": "mm",
    "press_msl": "hPa",
    "temp": "degC",
    "kin": "W m-2",
}

UNITS_MULT = {
    "precip": "*",
    "press_msl": "*",
    "temp": "+",
    "kin": "*",
}


class MyDumper(yaml.SafeDumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()


def check_directory(
    p: Path | str,
    root: Path | str = None,
):
    """_summary_."""
    p = Path(p)
    if not p.is_absolute():
        if root is None:
            raise ValueError(f"Path {p.as_posix()} is relative and root is None")
        p = Path(root, p)
    if not p.exists():
        p.mkdir(parents=True)
    return p


def posix_mount_check(
    target: Path,
    root: Path,
):
    """Check mount on unix systems."""
    target_mount = target.parts[1]
    root_mount = root.parts[1]
    if target_mount != root_mount:
        return False
    return True


def make_path_relative(
    target: Path | str,
    root: Path | str,
    max_depth: int = 5,
) -> Path:
    """Make a path relative to another.
    The target path will be made relative to the root path.
    This will be done when the following conditions are met:
    - `target` and `root` are on the same mount
    - directory depth difference does not exceed `max_depth`
    Parameters
    ----------
    target : Path | str
        The path that is to be made relative.
    root : Path | str
        The path that is used as the root againts which `target` is made relative.
    max_depth : int, optional
        The max difference allowed in depth between the `target` path and \
the `root` path, by default 5
    Returns
    -------
    Path
        Relative path (or absolute when conditions are not met).
    """
    # ensure pathlib objects
    target = Path(target)
    root = Path(root)
    # check for the mount
    same_mount = True
    if isinstance(target, WindowsPath):
        same_mount = target.drive == root.drive
    else:
        same_mount = posix_mount_check(target, root)
    if not same_mount:
        return target

    # get the relative path
    rel_path = Path(
        relpath(
            target,
            root,
        )
    )
    # If it exceeds the maximum depth
    if len(rel_path.parts) > max_depth:
        return target
    return rel_path