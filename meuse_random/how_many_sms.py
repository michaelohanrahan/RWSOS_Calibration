from glob import glob
from pathlib import Path

path = Path("/p/11209265-grade2023/wflow/RWSOS_Calibration/meuse_random/data/2-interim/calib_data")


sms = glob(str(path / "level1" / "**" / "staticmaps.nc"))

print(len(sms))
