plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 211, in _acquire_with_cache_info
    file = self._cache[self._key]
           ~~~~~~~~~~~^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\lru_cache.py", line 56, in __getitem__
    value = self._cache[key]
            ~~~~~~~~~~~^^^^^
KeyError: [<class 'netCDF4._netCDF4.Dataset'>, ('p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc',), 'r', (('clobber', True), ('diskless', False), ('format', 'NETCDF4'), ('persist', False)), '291c64af-4a50-4344-87bb-74685b854d48']

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<ipython-input-19-36031e1afa4c>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-19-36031e1afa4c>", line 66, in main
    obs_ds = xr.open_dataset(obs_data)
             ^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\api.py", line 588, in open_dataset
    backend_ds = backend.open_dataset(
                 ^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 645, in open_dataset
    store = NetCDF4DataStore.open(
            ^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 408, in open
    return cls(manager, group=group, mode=mode, lock=lock, autoclose=autoclose)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 355, in __init__
    self.format = self.ds.data_model
                  ^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 417, in ds
    return self._acquire()
           ^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\netCDF4_.py", line 411, in _acquire
    with self._manager.acquire_context(needs_lock) as root:
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 199, in acquire_context
    file, cached = self._acquire_with_cache_info(needs_lock)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\backends\file_manager.py", line 217, in _acquire_with_cache_info
    file = self._opener(*self._args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "src\netCDF4\_netCDF4.pyx", line 2464, in netCDF4._netCDF4.Dataset.__init__
  File "src\netCDF4\_netCDF4.pyx", line 2027, in netCDF4._netCDF4._ensure_nc_success
FileNotFoundError: [Errno 2] No such file or directory: 'p:\\p\\11209265-grade2023\\wflow\\RWSOS_Calibration\\meuse_random\\data\\1-external\\discharge_hourlyobs_smoothed.nc'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-2-325fc8c8eaa0>", line 221, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-2-325fc8c8eaa0>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-3-b952dce436d4>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-3-b952dce436d4>", line 75, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-3-b952dce436d4>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-3-b952dce436d4>", line 75, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-4-b8bc4c345491>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-4-b8bc4c345491>", line 76, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-4-b8bc4c345491>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-4-b8bc4c345491>", line 76, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-4-b8bc4c345491>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-4-b8bc4c345491>", line 76, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-2-4a687cb7aa89>", line 224, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-2-4a687cb7aa89>", line 75, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-3-37c55e9cd95a>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-3-37c55e9cd95a>", line 75, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-3-37c55e9cd95a>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-3-37c55e9cd95a>", line 75, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - cannot access local variable 'level_list' where it is not associated with a value
plot_final_model.py - ERROR - cannot access local variable 'level_list' where it is not associated with a value
plot_final_model.py - ERROR - cannot access local variable 'level_list' where it is not associated with a value
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-5-d22ea79c6f54>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-5-d22ea79c6f54>", line 63, in main
    level_list = [level.replace('level-1', 'base') for level in level_list]
                                                                ^^^^^^^^^^
UnboundLocalError: cannot access local variable 'level_list' where it is not associated with a value

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-5-d22ea79c6f54>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-5-d22ea79c6f54>", line 63, in main
    level_list = [level.replace('level-1', 'base') for level in level_list]
                                                                ^^^^^^^^^^
UnboundLocalError: cannot access local variable 'level_list' where it is not associated with a value

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-5-d22ea79c6f54>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-5-d22ea79c6f54>", line 63, in main
    level_list = [level.replace('level-1', 'base') for level in level_list]
                                                                ^^^^^^^^^^
UnboundLocalError: cannot access local variable 'level_list' where it is not associated with a value

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-6-4d549710d89f>", line 220, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-6-4d549710d89f>", line 72, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-6-4d549710d89f>", line 220, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-6-4d549710d89f>", line 72, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-6-4d549710d89f>", line 220, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-6-4d549710d89f>", line 72, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-6-4d549710d89f>", line 220, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-6-4d549710d89f>", line 72, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-7-b62ccfc0b985>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-7-b62ccfc0b985>", line 73, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-7-b62ccfc0b985>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-7-b62ccfc0b985>", line 73, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-7-b62ccfc0b985>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-7-b62ccfc0b985>", line 73, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-7-b62ccfc0b985>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-7-b62ccfc0b985>", line 73, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-7-b62ccfc0b985>", line 222, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-7-b62ccfc0b985>", line 73, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'
plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

plot_final_model.py - ERROR - Traceback (most recent call last):
  File "<ipython-input-8-a55ebd46f131>", line 223, in <module>
    ds = main(
         ^^^^^
  File "<ipython-input-8-a55ebd46f131>", line 74, in main
    da = xr.DataArray(
         ^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 471, in __init__
    coords, dims = _infer_coords_and_dims(data.shape, coords, dims)
                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 202, in _infer_coords_and_dims
    _check_coords_dims(shape, new_coords, dims_tuple)
  File "c:\Users\ohanrah\.conda\envs\hmt\Lib\site-packages\xarray\core\dataarray.py", line 132, in _check_coords_dims
    raise ValueError(
ValueError: conflicting sizes for dimension 'time': length 1 on the data but length 48 on coordinate 'time'

