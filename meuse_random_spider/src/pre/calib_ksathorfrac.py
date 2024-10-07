import copy
import os
import subprocess
from pathlib import Path

import tomli
import tomli_w

NEWLINE_CHAR = os.linesep

def main(
    ref_toml: Path | str,
    ref_shell: Path | str,
    factors: tuple | list,
):
    """_summary_."""
    p_toml = Path(ref_toml)
    p_shell = Path(ref_shell)

    with open(p_toml, "rb") as f:
        ref_toml = tomli.load(f)

    with open(p_shell, "r") as f:
        ref_shell = f.read().split("\n")

    for _f in factors:
        new_toml = copy.deepcopy(ref_toml)
        new_shell = copy.deepcopy(ref_shell)
        nname = f"ksat_{_f}"

        new_toml["dir_output"] = f"run_default/{nname}"
        new_toml["input"]["lateral"]["subsurface"]["ksathorfrac"] = {"value": _f}
        new_shell[1] = new_shell[1].replace("ref", nname)
        new_shell[2] = new_shell[2].replace("ref", nname)
        new_shell[-1] = new_shell[-1].replace("ref", nname)

        out_toml = p_toml.name.replace("ref", nname)
        out_shell = p_shell.name.replace("ref", nname)

        with open(Path(p_toml.parent, out_toml), "wb") as w:
            tomli_w.dump(new_toml, w)
        
        with open(Path(p_shell.parent, out_shell), "w") as w:
            w.write("\n".join(new_shell))
        
        subprocess.run(["dos2unix", Path(p_shell.parent, out_shell).as_posix()])


if __name__ == "__main__":
    factors = [5,10,20,25,50,75,100,200,250,500,750,1000]
    main(
        r"p:\1000365-002-wflow\tmp\usgs_wflow\models\MODELDATA_KING_500M\wflow_sbm_ref.toml",
        r"p:\1000365-002-wflow\tmp\usgs_wflow\jobs\CALIB\run_king_ref.sh",
        factors,
    )
    pass