import hydromt
import hydromt_wflow
from hydromt_wflow import WflowModel
from hydromt.log import logger
from helper import syscheck

drive = syscheck()

# Set the working directory
os.chdir(drive + r"\11209265-grade2023\wflow\RWSOS_Calibration")
logger.setLevel("INFO")
mod_root= drive + r"\11209265-grade2023\wflow\wflow_meuse_julia\wflow_meuse_202312"
new_root = 
# Initialize the model
model = WflowModel(root=mod_root, mode="build", name="meuse")

