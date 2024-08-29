import json
import os
import sys
from pathlib import Path


def false_index(
    oid: str | int,
    data: list,
):
    """_summary_"""
    _check = [
        oid in _d for _d in data
    ]
    if not any(_check):
        return None
    
    return _check.index(True) + 1


def sort_graph(
    p: Path | str,
):
    """_summary_"""
    with open(p, "r") as _r:
        data = json.load(_r)

    _count = 0
    _lvls = {}
    _lvl1 = [
        key for key, item in data.items() if len(item["_deps"]) == 0 
    ]
    for item in _lvl1:
        data.pop(item)
    _lvls[1] = _lvl1
    _rem = list(data.keys()) 
    
    while True:
        for item in _rem[:]:
            _check = [
                false_index(_dep, _lvls.values()) for _dep in data[item]["_deps"]
            ]
            if None in _check:
                continue
            _lvl = max(_check) + 1

            if _lvl not in _lvls:
                _lvls[_lvl] = []

            _lvls[_lvl].append(item)
            _rem.remove(item)
        
        if len(_rem) == 0:
            break

        if _count == 1000:
            sys.stdout.write(f"Could not solve the following id's: {_rem}\n")
            break
        _count += 1

    _out = {}
    for key, item in _lvls.items():
        dep = [f"level{key-1}"]
        if key-1 == 0:
            dep = []

        _out[f"level{key}"] = {
            "deps": dep,
            "elements": item,
        }
    
    return _out


def create_graph():
    pass


if __name__ == "__main__":
    p = Path("c:/CODING/NONPRODUCT/puget/res/king_graph.json")
    sort_graph(p)
    pass