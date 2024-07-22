#explanation of creating the health check

This is the precursor that creates all gauge location geojson
Filtering is done in QGIS and is noted in the Onenote 

the files are created in the meuse directory via:

CLI: 
'''
cd meuse
pixi shell
python src/pre/assess_discharge_data.py --freq D
'''

after QGIS sanitizing we save the interim files at
meuse\data\2-interim\QGIS\to_wflow\