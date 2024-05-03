import json
from pathlib import Path

import geopandas as gpd


def select_gauges(
    sel: Path | str,
    gauges: Path | str,
):
    """_summary_"""
    with open(sel, "r") as _r:
        data = json.load(_r)
    
    gauges = Path(gauges)
    out_dir = gauges.parent
    bname = gauges.stem
    suffix = gauges.suffix 

    ds = gpd.read_file(gauges)
    items = [
        int(item) for item in data.keys()
    ]
    ds = ds[ds.SiteNumber.isin(items)]

    ds.to_file(
        Path(out_dir, f"{bname}_calib{suffix}")
    )

    pass



if __name__ == "__main__":
    select_gauges(
        "c:/CODING/NONPRODUCT/puget/res/king_graph.json",
        "p:/1000365-002-wflow/tmp/usgs_wflow/models/GIS/gauges_king.gpkg",
    )