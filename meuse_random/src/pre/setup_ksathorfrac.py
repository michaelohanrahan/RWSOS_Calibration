import argparse
import os
from hydromt_wflow import WflowModel
from hydromt.data_catalog import DataCatalog
from hydromt.log import setuplog
import hydromt_wflow
import traceback
from icecream import ic
import shutil
import yaml

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
    logger.info(f'vars passed: {DRIVE}, {mod_root}, {mod_new_root}, {config_fn_in}, {config_fn_out}, {var}')
    mod = WflowModel(root=mod_root, 
                     data_libs=[f'{DRIVE}/wflow_global/hydromt_wflow/catalog.yml'],
                     config_fn=config_fn_in, 
                     mode='r', 
                     logger=logger)
    
    mod.read_config()
    mod.read_geoms()
    mod.read_grid()

    mod.set_root(mod_new_root, mode='w+')
    
    #the existing hydromt_wflow catalog meta is for windows, here it depends on the OS
    with open(f'{DRIVE}/wflow_global/hydromt_wflow/catalog.yml', 'r') as file:
        dc = yaml.safe_load(file)
    dc['meta']['root'] = f'{DRIVE}/wflow_global/hydromt_wflow'
    
    with open(f'{os.getcwd()}/catalog.yml', 'w') as file:
        yaml.dump(dc, file)
    dc = DataCatalog(f'{os.getcwd()}/catalog.yml')


    mod.setup_ksathorfrac(dc.get_rasterdataset('ksathorfrac')[f'{var}250'])
    key = [key for key in mod.grid.keys() if f'{var}250' in key][0]
    mod.set_grid(mod.grid.rename({key: f'ksathorfrac_{var}250'}))
    mod.config['input']['path_static'] = f'staticmaps/staticmaps.nc'
    mod.config['input']['lateral']['subsurface']['ksathorfrac'] = f'ksathorfrac_{var}250'
    mod.write_config(config_name=config_fn_out)
    mod.write_grid()
    if os.path.exists(os.path.join(mod_root, 'staticgeoms')):
        shutil.rmtree(os.path.join(mod_new_root, 'staticgeoms'))
        shutil.copytree(os.path.join(mod_root, 'staticgeoms'), os.path.join(mod_new_root, 'staticgeoms'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='WflowModel setup script')
    parser.add_argument('--DRIVE', required=True, help='Drive letter')
    parser.add_argument('--mod_root', required=True, help='Model root directory')
    parser.add_argument('--mod_new_root', required=True, help='New model root directory')
    parser.add_argument('--config_fn_in', required=True, help='Input configuration filename')
    parser.add_argument('--config_fn_out', required=True, help='Output configuration filename')
    parser.add_argument('--var', required=True, help='Variable prefix (e.g., BRT_ or RF_)')
    
    args = parser.parse_args()
    
    print(args.mod_root)
    try:
        main(args.DRIVE, 
            args.mod_root,
            args.mod_new_root,
            args.config_fn_in,
            args.config_fn_out,
            args.var)
    except Exception as e:
        print(f'Error: {e}')
        traceback.print_exc()