using Pkg
Pkg.activate(raw"P:\11209265-grade2023\wflow\wflow_meuse_julia\compare_fl1d_interreg\wflow_073")
Pkg.instantiate()
Pkg.status()

using Wflow

dirname = "run_fl1d_lakes"

config_str = string(dirname, ".toml")
config = Wflow.Config(config_str)
Wflow.run(config)
