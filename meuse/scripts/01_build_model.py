import hydromt
import hydromt_wflow
from hydromt_wflow import WflowModel
from hydromt.log import logger
import os 

os.chdir(r)

logger.setLevel("INFO")
mod_root= r"p:\11209265-grade2023\wflow\wflow_meuse_julia\wflow_meuse_202312"
new_root = 
# Initialize the model
model = WflowModel(root=mod_root, mode="build", name="meuse")

