import json
from itertools import product
from pathlib import Path
import pandas as pd
import numpy as np


# def create_set(
#     p: Path | str,
# ):
#     """
#     This function takes a path to a JSON file as input, reads the file, and creates a pandas DataFrame from the data.
#     The DataFrame is created by taking the Cartesian product of the "values" lists from the JSON data, and using the "short_name" fields as column names.
#     The function also returns a list of the keys in the JSON data, and a list of the "method" fields from the JSON data.

#     Parameters:
#     p (Path | str): The path to the JSON file to read.

#     Returns:
#     lnames (list): A list of the keys in the JSON data.
#     methods (list): A list of the "method" fields from the JSON data.
#     ds (DataFrame): A pandas DataFrame created from the JSON data.
#     """
#     with open(p, "r") as _r:
#         data = json.load(_r)

#     params = tuple(
#         product(*[item["values"] for item in data.values()]),
#     )
    
#     # Extract column names, 'nl, nf' is separated
#     columns = []
#     for item in data.values():
#         if ',' in item["short_name"]:
#             columns.extend(item["short_name"].split(','))
#         else:
#             columns.append(item["short_name"])
    
#     ds = pd.DataFrame(params, columns=columns[:-1])
#     ds[columns[-1]] = ds['nl']  # add 'nf' column with cell value the same as 'nl' column
    
#     lnames = list(data.keys())
#     last_element = lnames[-1]  # Separate the last element into two parts
#     split_elements = last_element.split(',')
#     lnames = lnames[:-1] + split_elements # Replace the last element with the new components
    
#     _methods = [
#         item["method"] for item in data.values()
#     ]
#     methods = _methods[:] + [_methods[-1]]  # add methods for N_Floodplein the same as N_Land
    
    
#     return lnames, methods, ds


def create_set(p: Path | str):
    """
    This function takes a path to a JSON file as input, reads the file, and creates a pandas DataFrame from the data.
    The DataFrame is created by taking the Cartesian product of the "values" lists from the JSON data, and using the "short_name" fields as column names.
    The function also returns a list of the keys in the JSON data, and a list of the "method" fields from the JSON data.

    Parameters:
    p (Path | str): The path to the JSON file to read.

    Returns:
    lnames (list): A list of the keys in the JSON data.
    methods (list): A list of the "method" fields from the JSON data.
    ds (DataFrame): A pandas DataFrame created from the JSON data.
    """
    with open(p, "r") as _r:
        data = json.load(_r)

    # Extract the cartesian product of all values
    params = list(product(*[item["values"] for item in data.values()]))

    # Extract column names, handling comma-separated short_names
    columns = []
    split_dict = {}
    for item in data.values():
        if ',' in item["short_name"]:
            split = item["short_name"].split(',')
            split_dict[split[0]] = split[1:]
            columns.append(split[0])
        else:
            columns.append(item["short_name"])
    
    # Create the DataFrame
    ds = pd.DataFrame(params, columns=columns)
   
    if len(split_dict) > 0:
        for col in ds.columns:
            if col in split_dict:
                for i in split_dict[col]:
                    ds[i] = ds[col]   

    # Extract the keys as lnames, handling comma-separated keys
    lnames = []
    for key in data.keys():
        lnames.extend(key.split(','))

    # Extract methods, handling cases where a method needs to be duplicated for split keys
    methods = []
    for item in data.values():
        method_count = len(item["short_name"].split(','))
        methods.extend([item["method"]] * method_count)

    return lnames, methods, ds


if __name__ == "__main__":
#     ds = create_set("c:/CODING/NONPRODUCT/puget/res/calib_recipe.json")
    fn = R"N:\My Documents\unix\documents\RWSoS\RWSOS_Calibration\meuse\config\MINIMAL_calib_recipe.json"
    lnames, methods, ds = create_set(fn)
    print(ds)