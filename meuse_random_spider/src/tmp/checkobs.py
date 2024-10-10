import xarray as xr 
from pathlib import Path
import os 

import pandas as pd
import ast

# Sample DataFrame
data = {
    'gauges': [1, 2, 3],
    'Top_1': ["{'param1': 10, 'param2': 20}", "{'param1': 30, 'param2': 40}", "{'param1': 50, 'param2': 60}"]
}
params_ds = pd.DataFrame(data).set_index('gauges')

# Print the DataFrame
print("Original DataFrame:")
print(params_ds)

# Access the string in the DataFrame and convert it to a dictionary
gauge = 1
param_string = params_ds.loc[gauge, "Top_1"]
param_set = ast.literal_eval(param_string)

# Print the converted dictionary
print("\nConverted Dictionary:")
print(param_set)

# Example usage in a loop
for key, value in param_set.items():
    print(f"Key: {key}, Value: {value}")