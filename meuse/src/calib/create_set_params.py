import json
from itertools import product
from pathlib import Path

import pandas as pd


def create_set(
    p: Path | str,
):
    """_summary_"""
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


if __name__ == "__main__":
    ds = create_set("c:/CODING/NONPRODUCT/puget/res/calib_recipe.json")
    # ds = create_set(r'p:\11209265-grade2023\wflow\RWSOS_Calibration\meuse\config\calib_recipe.json')