import xarray as xr 

sm_old = xr.open_dataset("staticmaps/staticmaps_routing_cal_11.nc")

sm_int = xr.open_dataset(r"P:\11209265-grade2023\wflow\wflow_meuse_julia\compare_fl1d_interreg\interreg\staticmaps\staticmaps.nc")

add_vars = ['wflow_lakeareas', 
            'wflow_lakelocs',
            'LakeArea', 
            'LakeAvgLevel', 
            'LakeAvgOut',
            'Lake_b',
            'Lake_e',
            'LakeStorFunc',
            'LakeOutflowFunc',
            'LakeThreshold',
            'LinkedLakeLocs'
            ]

for var in add_vars:
    sm_old[var] = sm_int[var]

sm_old.to_netcdf("staticmaps/staticmaps_add_lakes.nc")