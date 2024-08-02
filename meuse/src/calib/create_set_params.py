import json
from itertools import product
from pathlib import Path
import pandas as pd


def create_set(
    p: Path | str,
):
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

    params = tuple(
        product(*[item["values"] for item in data.values()]),
    )
    
    # Extract column names, 'nl, nf' is separated
    columns = []
    for item in data.values():
        if ',' in item["short_name"]:
            columns.extend(item["short_name"].split(','))
            SPLIT = True
        else:
            columns.append(item["short_name"])
            SPLIT = False
    # print(columns)
    # print(params)
    if SPLIT:
        ds = pd.DataFrame(params, columns=columns[:-1])
        ds[columns[-1]] = ds['nl']
        lnames = list(data.keys())
        last_element = lnames[-1]  # Separate the last element into two parts
        split_elements = last_element.split(',')
        lnames = lnames[:-1] + split_elements # Replace the last element with the new components
        
        _methods = [
            item["method"] for item in data.values()
        ]
        methods = _methods[:] + [_methods[-1]]  # add methods for N_Floodplein the same as N_Land
    else:
        ds = pd.DataFrame(params, columns=columns)
        lnames = list(data.keys())
        methods = [
            item["method"] for item in data.values()
        ]

    # check if lnames, methods and ds.columns have the same length
    # if len(lnames) == len(methods) == len(ds.columns.values):
    #     print("Perfect!")
    # else:
    #     print("Error: lnames, methods and ds.columns do not have the same length!")
    
    return lnames, methods, ds


if __name__ == "__main__":
#     ds = create_set("c:/CODING/NONPRODUCT/puget/res/calib_recipe.json")
    fn = R"c:\Users\deng_jg\work\05wflowRWS\RWSOS_Calibration\meuse\config\calib_recipe.json"
    lnames, methods, ds = create_set(fn)
    print(ds)