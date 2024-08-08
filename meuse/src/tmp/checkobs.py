import xarray as xr 
from pathlib import Path
import os 

import pandas as pd
import ast

d= {   
    "ksathorfrac_BRT_250": { 
        "short_name": "ksat", 
        "values": [0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0],
        "method": "mult"
    },
    "f_":{
        "short_name": "f",
        "values": [0.5, 0.75, 1.0, 1.25, 1.5],
        "method": "mult"
    },
    "RootingDepth_obs_15": {
        "short_name": "rd",
        "values": [0.25, 0.5, 1, 1.5, 2.0, 2.5],
        "method": "mult"
    },
    "SoilThickness_manual_cal": {
        "short_name": "st",
        "values": [0.5, 0.75, 1.0, 1.25, 1.5],
        "method": "mult"
    },
    "N_River": {
        "short_name": "nr",
        "values": [0.5, 0.75, 1.0, 1.25, 1.5],
        "method": "mult"
    }, 
    "MaxLeakage_manual_cal": {
        "short_name": "ml",
        "values": [0, 0.2, 0.6],
        "method": "add"
    }, 
    "N,N_Floodplain": {
        "short_name": "nl,nf",
        "values": [0.8, 1.0, 1.2],
        "method": "mult"
    }
}
vals = 1

for key, r in d.items():
    vals *= len(r["values"])
print(vals)