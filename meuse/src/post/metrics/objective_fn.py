import numpy as np


# ------------ Calculate NSE and NSE log --------------
def calculate_nse_and_log_nse(observed, modelled):
    """
    Calculates the Nash-Sutcliffe Efficiency (NSE) and Log-Nash-Sutcliffe Efficiency (NSE_log) 
    between observed and modelled data.
    
    Parameters:
    observed (array-like): Array of observed data.
    modelled (array-like): Array of modelled data.
    
    Returns:
    nse (float): Nash-Sutcliffe Efficiency.
    nse_log (float): Log-Nash-Sutcliffe Efficiency.
    """
    observed = np.array(observed)
    modelled = np.array(modelled)

    obs_mean = np.mean(observed)

    numerator = np.sum((observed - modelled) ** 2)

    denominator = np.sum((observed - obs_mean) ** 2)

    nse = 1 - (numerator / denominator)

    log_numerator = np.sum((np.log(observed + 1) - np.log(modelled + 1)) ** 2)
    log_denominator = np.sum((np.log(observed + 1) - np.log(obs_mean + 1)) ** 2)
    nse_log = 1 - (log_numerator / log_denominator)

    return nse, nse_log

