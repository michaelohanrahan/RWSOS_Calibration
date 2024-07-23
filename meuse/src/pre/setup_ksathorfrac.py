import argparse
import os
from hydromt_wflow import WflowModel
from hydromt.data_catalog import DataCatalog
from hydromt.log import setuplog
import hydromt_wflow

def main(DRIVE:str, 
         mod_root:str, 
         mod_new_root:str, 
         config_fn_in:str, 
         config_fn_out:str, 
         var:str):
    """
    SUMMARY
    """
    logger = setuplog("wflow", log_level=10)
    logger.info(f'Hydromt version: {hydromt_wflow.__version__}')

    mod = WflowModel(root=mod_root, 
                     config_fn=config_fn_in, 
                     mode='r', 
                     logger=logger)
    
    mod.read_config()
    mod.read_geoms()
    mod.read_grid()

    mod.set_root(mod_new_root, mode='w+')

    dc = DataCatalog()
    dc = dc.from_yml(urlpath=f'{DRIVE}/wflow_global/hydromt_wflow/catalog.yml')

    # Problematic ksatname
    rename = {f'{DRIVE}/wflow_global/hydromt_wflow/soil/ksathorfrac.zarr_{var}250': f'ksathorfrac_{var}250'}

    mod.setup_ksathorfrac(dc['ksathorfrac'].path, variable=f'{var}250')

    # Rename the grid variable in the dataset
    old_name = list(rename.keys())[0]
    new_name = list(rename.values())[0]

    mod.set_grid(mod.grid.rename_vars({old_name: new_name}))

    mod.config['input']['path_static'] = f'staticmaps/staticmaps.nc'
    mod.config['input']['lateral']['subsurface']['ksathorfrac'] = new_name
    mod.write_config(config_name=config_fn_out)
    mod.write_grid()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='WflowModel setup script')
    parser.add_argument('--DRIVE', required=True, help='Drive letter')
    parser.add_argument('--mod_root', required=True, help='Model root directory')
    parser.add_argument('--mod_new_root', required=True, help='New model root directory')
    parser.add_argument('--config_fn_in', required=True, help='Input configuration filename')
    parser.add_argument('--config_fn_out', required=True, help='Output configuration filename')
    parser.add_argument('--var', required=True, help='Variable prefix (e.g., BRT_ or RF_)')

    args = parser.parse_args()
    main(args.DRIVE, 
         args.mod_root,
         args.mod_new_root,
         args.config_fn_in,
         args.config_fn_out,
         args.var)