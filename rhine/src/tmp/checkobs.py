import xarray as xr 
from pathlib import Path
import os 
loc = Path("p:/","11209265-grade2023/wflow/RWSOS_Calibration","rhine","data/3-input","inmaps/FINAL_forcing-genRE_hourly_*.nc").as_posix()
print(loc)
#list files in dir 
files = os.listdir(os.path.dirname(loc))
print(files)
ds = xr.open_mfdataset(loc)
ds
