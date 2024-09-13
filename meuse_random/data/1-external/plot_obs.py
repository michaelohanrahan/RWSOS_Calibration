import xarray as xr
import matplotlib.pyplot as plt
import os 

os.chdir(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse_random\data\1-external')
ds = xr.open_dataset(os.path.join(os.getcwd(),'discharge_hourlyobs_smoothed.nc'))
t0='2017-06-25'
t1='2017-07-30'
fig,ax = plt.subplots(1,1, figsize=(10,10))
ds.sel(runs='Obs.', wflow_id=16, time=slice(t0,t1)).Q.plot(ax=ax, c='orange', alpha=0.5, label='obs')
ds.sel(runs='HBV', wflow_id=16, time=slice(t0,t1)).Q.plot(ax=ax, c='blue', alpha=0.5, label='hbv')
plt.legend()
plt.show()