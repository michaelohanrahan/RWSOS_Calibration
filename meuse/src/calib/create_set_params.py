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
    
    ds = pd.DataFrame(
        params,
        columns=[
            item["short_name"] for item in data.values()
        ]
    )

    lnames = list(data.keys())
    methods = [
        item["method"] for item in data.values()
    ]
    return lnames, methods, ds


# if __name__ == "__main__":
# #     ds = create_set("c:/CODING/NONPRODUCT/puget/res/calib_recipe.json")
#     ds = create_set(r'c:\git\RWSOS_Calibration\meuse\config\calib_recipe.json')
#     print(ds)