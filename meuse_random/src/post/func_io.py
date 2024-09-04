import pandas as pd
from datetime import datetime

def read_filename_txt(filename):
    dfm = pd.read_csv(filename, sep=';', parse_dates=[0], header=None,skiprows=36,na_values=-999, encoding_errors='ignore')
    dfm=dfm.drop([1], axis=1)
    dfm=dfm.drop([0], axis=0)
    dfm.index=dfm[0]
    dfm.index.name=None
    dfm=dfm.drop([0], axis=1)
    return dfm


def read_lakefile(filename):
    dateparse = lambda x: datetime.strptime(x[17:34], '%Y.%m.%d %M:%S')
    dfm = pd.read_csv(filename, sep=';', parse_dates=[1], header=None,skiprows=7,na_values=-999,date_parser=dateparse)
    dfm=dfm.drop([0], axis=1)
    dfm.index=dfm[1]
    dfm=dfm.drop([1], axis=1)
    dfm=dfm.resample('D', label='right').mean()
    return dfm