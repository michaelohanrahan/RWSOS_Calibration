import os
import pandas as pd
import numpy as np
import xarray as xr
import plotly.graph_objects as go
from datetime import datetime
# from hydromt_wflow import WflowModel
from metrics.objective_fn import calculate_nse_and_log_nse
import traceback
import matplotlib.pyplot as plt
from datetime import datetime
from icecream import ic


def plot_peaks_ts(ds:xr.Dataset,  
            df_GaugeToPlot:pd.DataFrame, 
            start:datetime, 
            end:datetime, 
            Folder_plots:str, 
            color_dict:dict,
            run_keys:list=None,
            peak_dict:dict=None,
            savefig:bool=False,
            font:dict={'family': 'serif', 'size': 16},
            translate:dict=None,
            id_key:str='wflow_id'
            )->None:
    
    '''
    ds: xarray dataset that contains modeled results for all runs and observation data
        requires that the observations are indexed as 'Obs.'
    '''
    
    if run_keys is None or len(run_keys) == 0:
        run_keys = ds.runs.values
    
    for station_id in ds[id_key].values:
        # try:
        station_name = df_GaugeToPlot.loc[df_GaugeToPlot[id_key]==station_id, 'location'].values[0]

        print('Attempt html peak timing plotting for (station_name, station_id): ', (station_name, station_id))

        # try:
        fig = go.Figure()
        
        for run in ds.runs.values:
            
            obs = ds.sel(time=slice(start, end), runs='Obs.', wflow_id=station_id)
            obs = obs.ffill('time') #observations are prone to nans
            sim = ds.sel(time=slice(start, end), runs=run, wflow_id=station_id)
            
            # print('len obs, len sim', len(obs.Q), len(sim.Q))
        
            # handle nan values in the 
            mask = ~obs.Q.isnull() & ~sim.Q.isnull()
            
            obs_filtered = obs.Q[mask]
            sim_filtered = sim.Q[mask]
            
            
            if len(obs_filtered) == len(sim_filtered):
                nse, nse_log = calculate_nse_and_log_nse(obs_filtered, sim_filtered)
                # print('nse, nse_log', (nse, nse_log))
                
            else:
                nse, nse_log = np.nan, np.nan
                # print('nse, nse_log = np.nan, np.nan')    
                # print(f'len(obs_filtered) {len(obs_filtered)} != len(sim_filtered) {len(sim_filtered)}')
            
            if run == 'Obs.':
                label = f'{run}'
                
                obs_peaks = peak_dict[station_id][run]['peaks']
                
                fig.add_trace(go.Scatter(
                    x=obs.time.values,
                    y=obs.Q.values,
                    mode='lines',
                    name=label,
                    line=dict(color=color_dict[str(run)])
                ))
                
                
                fig.add_trace(go.Scatter(
                    x=obs.time.sel(time=obs_peaks).values,
                    y=obs.Q.sel(time=obs_peaks).values,
                    mode='markers',
                    name='Obs. Peaks',
                    marker=dict(color=color_dict[str(run)], size=8)
                ))
                        
            else:
                try:
                    label = f"{run}: {np.mean(peak_dict[station_id][run]['timing_errors']):.2f} +/- {np.std(peak_dict[station_id][run]['timing_errors']):.2f}h" + "\nNSE:" + f"{nse:.3f}" + " NSE<sub>log</sub>:" + f"{nse_log:.3f}"
                    
                    fig.add_trace(go.Scatter(
                        x=sim.time.values,
                        y=sim.Q.values,
                        mode='lines',
                        name=label,
                        line=dict(color=color_dict[str(run)])
                    ))
                    
                    sim_peaks = peak_dict[station_id][run]['peaks']

                    fig.add_trace(go.Scatter(
                        x=sim.time.sel(time=sim_peaks, method='nearest').values,
                        y=sim.Q.sel(time=sim_peaks, method='nearest').values,
                        mode='markers',
                        name='Sim. Peaks',
                        marker=dict(color=color_dict[str(run)], size=8)
                    ))
                except Exception as e:
                    print('Error in plotting peak timing for run, station:', run)
                    print(e)
                    traceback.print_exc()
                    None
                    
        fig.update_layout(
            title=f'{station_name} (id: {station_id})',
            xaxis_title='Date (hourly timestep)',
            yaxis_title='Discharge (m<sup>3</sup>s<sup>-1</sup>)',
            font=font
        )

        if savefig == True:
            interactive_folder = os.path.join(Folder_plots, 'interactive')
            os.makedirs(interactive_folder, exist_ok=True)
            filename = f'PeakTiming_{station_name}_{station_id}_{start.year}{start.month}{start.day}-{end.year}{end.month}{end.day}.html'
            filepath = os.path.join(interactive_folder, filename)
            fig.write_html(filepath)
            if os.path.exists(filepath):
                print(f'Peak timing plot for {station_name} (id: {station_id}) saved as\n{filename}')
            else:
                print(f'Failed to save {filepath}')
        else:
            fig.show()
        
    return None

#======================= Compute peak timing errors for all runs (and plot analysis results) =============
def peak_timing_for_runs(ds:xr.Dataset,
                         df_GaugeToPlot:pd.DataFrame,
                         folder_plots:str, 
                         peak_dict:dict, 
                         plotfig:bool=False, 
                         savefig:bool=False) -> None:
    '''
    Plotting peak timing errors for all runs and stations in the dataset, displayed as distributions of timing errors and peak errors.
    ds: xarray dataset that contains modeled results for all runs and observation data
    df_GaugeToPlot: pandas dataframe that contains the information about the wflow_id to plot,
    
    '''
    
    for i, station_id in enumerate(ds.wflow_id.values):
        # try:
        station_name = df_GaugeToPlot.loc[df_GaugeToPlot['wflow_id']==station_id, 'location'].values[0]
        
        # try:
        # select a station
        ds_sub = ds.sel(wflow_id=station_id)

        # get obs data
        obs = ds_sub.sel(runs='Obs.').Q
        
        # plot figures
        if plotfig == True:
            
            peak_dict_sid = peak_dict[station_id]
            keys = list(peak_dict_sid.keys())
            markers = ['o', 's', '^', 'D', 'x', '*', '+', 'v', '<', '>']
            
            fig = plt.figure(figsize=(15, 10))  # Wider figure to accommodate side-by-side bar plots
            
            # Scatter plot of Qobs vs timing error
            ax1 = plt.subplot2grid((2, 2), (0, 0))  # Scatter plot spans the first row

            for (run, data), marker in zip(peak_dict_sid.items(), markers):
                if run != 'Obs.':
                    peaks = data['peaks']
                    timing_errors = data['timing_errors']
                    ax1.scatter(obs.sel(time=peaks), timing_errors, marker=marker, label=run)
                    # else:
                    #     print(f"Data contains non-finite values: {data}")

            ax1.legend()
            ax1.set_xlabel('Observed Q ($m s^{-1}$)')
            ax1.axhline(0, color='black', alpha=0.5, linestyle='--')
            ax1.set_ylabel('Timing Error ($h$)')
            ax1.set_title(f'Meuse at {station_name} (id: {station_id}) - Scatter plot of Qobs vs timing error')

            ax2 = plt.subplot2grid((2, 2), (1, 0))  # First plot in second row

            colors = ['skyblue' if key != 'Obs.' else 'grey' for key in keys]

            means = []

            # Set x-axis labels
            ax2.set_xticklabels(keys, rotation=15)

            for i, (run, data) in enumerate(peak_dict_sid.items()):

                data = data['timing_errors']
                
                if run == 'Obs.':
                    mean = 0
                    means.append(mean)
                    data = np.zeros(10)
                else:
                    mean = np.nanmean(data, axis=0)
                
                means.append(mean)
                
                # Create a boxplot
                bplot = ax2.boxplot(data, positions=[i+1], patch_artist=True, notch=True, vert=1, showfliers=False)  # Set showfliers=False to remove outliers

            # Draw a line between the means
            ax2.plot(range(1, len(means) + 1), means, color='r', linestyle='--')
            ax2.grid(axis='y')

            # Set colors for each box
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)

            [ax2.text(i+1, mean, round(mean, 1), ha='center', va='bottom', color='white', fontsize=13) for i, mean in enumerate(means)]
            [ax2.text(i+1, mean, round(mean, 1), ha='center', va='bottom', fontsize=11) for i, mean in enumerate(means)]

            # Set x-axis labels
            # ax2.set_xticklabels(keys, rotation=15)
                
            ax2.set_xlabel('Run')
            ax2.set_ylabel('Peak Timing (h)')
            ax2.set_title('Peak Timing Per Run')
            ax2.axhline(0, color='red', alpha=0.3, linestyle='--')

            # Mean absolute percentage peak error bar plot
            ax3 = plt.subplot2grid((2, 2), (1, 1))  # Second plot in second row
            peak_mapes = [peak_dict_sid[key]['peak_mape'] for key in keys]
            bars = ax3.bar(keys, peak_mapes, color=colors)  # add yerr parameter

            
            # Add the data values on top of the bars
            for bar in bars:
                yval = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2, 
                        yval + 1,  # Add the error and a small offset from the top of the bar
                        round(yval, 2), 
                        ha='center', 
                        va='bottom')
            
            # ax3.set_ylim(0, max(peak_mapes) + 5)
            ax3.set_ylim(0, np.nanmax(peak_mapes)*1.2)
            ax3.set_xticklabels(list(peak_dict_sid.keys()), rotation=15)
            ax3.set_xlabel('Run')
            ax3.set_ylabel('MAPE (100%)')
            ax3.set_title('Mean Absolute Percentage Peak error (MAPE, Instaneous Discharge)')

            # Adjust layout to prevent overlap
            plt.tight_layout()
            # plt.show()
            
            #define start and end as the first and last timestamp in the ds
            start = pd.to_datetime(ds.time.min().values.astype(datetime))
            end = pd.to_datetime(ds.time.max().values.astype(datetime))
            
            
            if savefig:
                timeseries_folder = os.path.join(folder_plots, 'Event_Timing_Metrics')
                os.makedirs(timeseries_folder, exist_ok=True)
                filename = f'PeakTimingMetrics_{station_name}_{station_id}_{start.year}{start.month}{start.day}-{end.year}{end.month}{end.day}.png'
                plt.savefig(os.path.join(timeseries_folder, filename), dpi=300)
            
            # if i == 0:
            #     break

        if savefig == True and plotfig == False:
            print('plotfig is False, no figure saved.')
        else:
            None
        
            
        # except Exception as e:
        #     print(f'An error occurred for {station_name} (id: {station_id}): {str(e)}')
        #     traceback.print_exc()



#======================= Peak timing distribution =======================

def plot_peak_timing_distribution(run_keys,
                                  peak_dict,
                                  color_dict,
                                  Folder_plots):

    fig, axs = plt.subplots(5, 1, figsize=(10, 24), sharex=True)  # Create 4 subplots

    bin_edges = np.arange(-70, 90, 10)

    # Plot combined histogram and KDE in the top subplot
    for var in run_keys:
        if var != 'Obs.':
            sns.ecdfplot(peak_dict[16][var]['timing_errors'], 
                        ax=axs[0], label=var, color=color_dict[var],
                        stat='percent', linewidth=2.5)

    axs[0].set_title('Empirical Cumulative Distribution of Relative Timing Data', fontsize=16)
    axs[0].set_ylabel('Percent (%)', fontsize=14)
    axs[0].grid(True, linestyle='--', alpha=0.6)
    axs[0].legend()

    # Plot separate histograms in the subsequent subplots
    for i, var in enumerate([v for v in run_keys if v != 'Obs.'], start=1):
        # if i == 4:
        #     break
        hist, _ = np.histogram(peak_dict[16][var]['timing_errors'], bins=bin_edges)
        axs[i].hist(bin_edges[:-1], bin_edges, weights=hist, edgecolor='black', linewidth=1.2, label=var, alpha=0.5, color=color_dict[var])
        axs[i].set_title(f'Histogram of {var}', fontsize=16)
        axs[i].set_xlabel('Lead <-- Value --> Lag', fontsize=14)
        axs[i].set_ylabel('Frequency', fontsize=14)
        axs[i].grid(True, linestyle='--', alpha=0.6)
        axs[i].set_ylim([0, 65])
        

    plt.tight_layout()
    plt.savefig(os.path.join(Folder_plots, 'Peak_timing_distribution.png'), dpi=600)
