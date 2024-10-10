import argparse
import os
from hydromt_wflow import WflowModel
from hydromt.data_catalog import DataCatalog
from hydromt.log import setuplog
import hydromt_wflow
import traceback

dc = DataCatalog('p:/wflow_global/hydromt_wflow/catalog.yml')
dc.get_rasterdataset('ksathorfrac')['BRT_250'].name


