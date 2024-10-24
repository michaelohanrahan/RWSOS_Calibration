import os
from pathlib import Path
from datetime import datetime

import hydromt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import xarray as xr
import pandas as pd
from tqdm import tqdm
import time


def rsquared(x, y):
    """Return R^2 where x and y are array-like, handling NaN values."""
    # Remove NaN values
    mask = ~np.isnan(x) & ~np.isnan(y)
    x_clean = x[mask]
    y_clean = y[mask]

    # Check if we have enough data points after removing NaNs
    if len(x_clean) < 2:
        return np.nan  # Not enough data to calculate R^2

    # Calculate R^2
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x_clean, y_clean)
    return r_value**2

def mm7q(
    ds: xr.Dataset | xr.DataArray,
    dry_month: list,
    window: int,
):  
    # calculate the MM7Q for each month
    _mm7q = ds.rolling(time=window).mean().resample(time='M').min('time').compute()
    
    # select out the MM7Q for selected dry months
    dry_month_start = dry_month[0]
    dry_month_end = dry_month[-1]
    months = _mm7q['time'].dt.month
    mm7q_dry_month = _mm7q.sel(time=_mm7q['time'].where((months>=dry_month_start)
                                                        &(months<=dry_month_end),
                                                        drop=True))
    
    return mm7q_dry_month


def process_labels_and_colors(dsq, labels, colors):
    """
    Processes the labels and colors for plotting based on the input dataset.

    This function checks if the labels are provided as a dictionary or a list. If a dictionary, it maps the old labels to new labels and updates the dataset and colors accordingly. If a list, it filters out labels that are not present in the dataset. It also ensures 'Obs.' is included if present in the dataset.

    Parameters:
    - dsq (xr.Dataset): The dataset containing the data to be plotted.
    - labels (dict or list): The labels to be used for plotting. If a dictionary, maps old labels to new labels. If a list, filters out labels not present in the dataset.
    - colors (dict or list): The colors to be used for plotting. If a dictionary, maps old labels to colors. If a list, filters out colors not corresponding to labels present in the dataset.

    Returns:
    - dsq (xr.Dataset): The updated dataset with new labels if necessary.
    - labels (list): The list of labels to be used for plotting.
    - colors (list): The list of colors to be used for plotting.
    """
    if isinstance(labels, dict):
        available_labels = [label for label in labels.keys() if label in dsq.runs.values]
        colors_new = {}
        new_dsq = xr.Dataset()
        for old_label, new_label in labels.items():
            if old_label in dsq.runs.values:
                new_dsq[new_label] = dsq['Q'].sel(runs=old_label)
                colors_new[new_label] = colors[old_label]
        if 'Obs.' in dsq.runs.values:
            new_dsq['Obs.'] = dsq['Q'].sel(runs='Obs.')
        colors = colors_new
        new_dsq = new_dsq.to_array(dim='runs').to_dataset(name='Q')
        dsq = new_dsq
        labels = list(labels.values())
        colors = [colors[label] for label in labels]
        print("Plotting for the following labels: ", labels)
        print("With colors: ", colors)
        print(dsq)
    else:
        labels = [label for label in labels if label in dsq.runs.values]
        print("Plotting for the following labels: ", labels)
    return dsq, labels, colors

def plot_signatures(
    dsq,
    labels: list[str] | dict[str, str],
    colors,
    Folder_out,
    station_name,
    station_id,
    lw=0.8,
    fs=8,
    save=False,
    legend_threshold=15,
):
    if xr.infer_freq(dsq.time).lower() == "D":
        window = 7
    elif xr.infer_freq(dsq.time).lower() == "h":
        window = 7 * 24
    elif xr.infer_freq(dsq.time).lower() == "3h":
        window = 7 * int(24 / 3)

    dsq, labels, colors = process_labels_and_colors(dsq, labels, colors)
    # first calc some signatures
    dsq["metrics"] = ["KGE", "NSE", "NSElog", "NSElog_mm7q", "RMSE", "MSE"]
    dsq["performance"] = (
        ("runs", "metrics"),
        np.zeros((len(dsq.runs), len(dsq.metrics))) * np.nan,
    )

    _tmp = dsq.sel(time=dsq.time[~dsq.Q.sel(runs="Obs.").isnull()].values)
    if len(_tmp.time) == 0:
        skip_obs = True

    else:
        skip_obs = False
        dsq = dsq.sel(time=dsq.time[~dsq.Q.sel(runs="Obs.").isnull()].values)
    # Create a dataframe to store the metrics
    metrics_df = pd.DataFrame(columns=['Station', 'Run', 'KGE', 'NSE', 'NSElog', 'NSElog_mm7q', 'RMSE', 'MSE', 'R2', 'R2_max_annual', 'R2_nm7q'])

    # perf metrics for single station
    if not skip_obs:
        for label in labels:
            # nse
            print("Run: ", label)
            nse = hydromt.stats.nashsutcliffe(
                dsq["Q"].sel(runs=label), dsq["Q"].sel(runs="Obs.")
            )
            dsq["performance"].loc[dict(runs=label, metrics="NSE")] = nse
            # nse logq
            nselog = hydromt.stats.lognashsutcliffe(
                dsq["Q"].sel(runs=label), dsq["Q"].sel(runs="Obs.")
            )
            dsq["performance"].loc[dict(runs=label, metrics="NSElog")] = nselog
            # kge
            kge = hydromt.stats.kge(dsq["Q"].sel(runs=label), dsq["Q"].sel(runs="Obs."))
            dsq["performance"].loc[dict(runs=label, metrics="KGE")] = kge["kge"]
            # rmse
            rmse = hydromt.stats.rmse(
                dsq["Q"].sel(runs=label), dsq["Q"].sel(runs="Obs.")
            )
            dsq["performance"].loc[dict(runs=label, metrics="RMSE")] = rmse
            # mse
            mse = hydromt.stats.mse(dsq["Q"].sel(runs=label), dsq["Q"].sel(runs="Obs."))
            dsq["performance"].loc[dict(runs=label, metrics="MSE")] = mse
            # nselog_mm7q
            sim_mm7q = mm7q(dsq["Q"].sel(runs=label), dry_month=[6, 7, 8, 9, 10], window=window)
            obs_mm7q = mm7q(dsq["Q"].sel(runs="Obs."), dry_month=[6, 7, 8, 9, 10], window=window)
            NSElog_mm7q = hydromt.stats.lognashsutcliffe(sim_mm7q, obs_mm7q)
            dsq["performance"].loc[dict(runs=label, metrics="NSElog_mm7q")] = NSElog_mm7q
            print(f"Run: {label}, NSE: {nse.values}, NSElog: {nselog.values}, KGE: {kge['kge'].values}, RMSE: {rmse.values}, MSE: {mse.values}, NSElog_mm7q: {NSElog_mm7q.values}") 

            # Calculate R2 scores
            r2_score = rsquared(dsq["Q"].sel(runs="Obs."), dsq["Q"].sel(runs=label))
            r2_score_max_annual = rsquared(dsq_max["Q"].sel(runs="Obs."), dsq_max["Q"].sel(runs=label))
            r2_score_nm7q = rsquared(dsq_nm7q["Q"].sel(runs="Obs."), dsq_nm7q["Q"].sel(runs=label))

            # Add metrics to the dataframe
            metrics_df = metrics_df.append({
                'Station': station_name,
                'Run': label,
                'KGE': kge['kge'].values[0],
                'NSE': nse.values,
                'NSElog': nselog.values,
                'NSElog_mm7q': NSElog_mm7q.values,
                'RMSE': rmse.values,
                'MSE': mse.values,
                'R2': r2_score,
                'R2_max_annual': r2_score_max_annual,
                'R2_nm7q': r2_score_nm7q
            }, ignore_index=True) 
    # needed later for sns boxplot
    #     df_perf = pd.DataFrame()
    #     for label in [label_00, label_01]:
    #         df = dsq['performance'].sel(runs = label, metrics = ['NSE', 'NSElog', 'KGE']).to_dataframe()
    #         df_perf = pd.concat([df,df_perf])

    fig = plt.figure(
        "signatures", clear=True, figsize=(24 / 2.54, 38 / 2.54), tight_layout=True
    )
    # fig, axes = plt.subplots(5, 2, figsize=(16 / 2.54, 22 / 2.54))
    axes = fig.subplots(5, 2)
    axes = axes.flatten()

    fig.suptitle(station_name)

    # daily against each other axes[0]
    print("Labels: ", labels)
    print("Colors: ", colors)
    if not skip_obs:
        for label, color in zip(labels, colors):
            print("Plotting: ", label)
            print("Color: ", color)
            axes[0].plot(
                dsq["Q"].sel(runs="Obs."),
                dsq["Q"].sel(runs=label),
                marker="o",
                linestyle="None",
                linewidth=lw,
                label=label,
                color=color,
                markersize=3,
            )

    max_y = np.round(dsq["Q"].max().values)
    axes[0].plot([0, max_y], [0, max_y], color="0.5", linestyle="--", linewidth=1)
    axes[0].set_xlim([0, max_y])
    axes[0].set_ylim([0, max_y])
    axes[0].set_ylabel("Simulated Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[0].set_xlabel("Observed Q (m$^3$s$^{-1}$)", fontsize=fs)
    #     axes[0].legend(frameon=True, fontsize = fs, )
    if len(labels) + 1 < legend_threshold:
        axes[0].legend(
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc="lower left",
            ncol=2,
            mode="expand",
            borderaxespad=0.0,
            fontsize=fs,
        )

    if len(labels) + 1 < legend_threshold and not skip_obs:
        text = ""
        for label in labels:
            r2_score = rsquared(dsq["Q"].sel(runs="Obs."), dsq["Q"].sel(runs=label))
            text += f"R$_2$ {label} = {r2_score:.2f} \n"
        axes[0].text(
            0.01,
            0.99,
            text,
            fontsize=fs,
            ha="left",
            va="top",
            transform=axes[0].transAxes,
        )

    # streamflow regime axes[1]
    for label, color in zip(labels, colors):
        dsq["Q"].sel(runs=label).groupby("time.month").mean("time").plot(
            ax=axes[1], linewidth=lw, label=label, color=color
        )

    if not skip_obs:
        dsq["Q"].sel(runs="Obs.").groupby("time.month").mean("time").plot(
            ax=axes[1], linewidth=lw, label="Obs.", color="k", linestyle="--"
        )
    axes[1].tick_params(axis="both", labelsize=fs)
    axes[1].set_ylabel("Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[1].set_xlabel("month", fontsize=fs)
    axes[1].set_title("")
    axes[1].set_xticks(np.arange(1, 13))
    if len(labels) + 1 < legend_threshold:
        axes[1].legend(
            bbox_to_anchor=(0.0, 1.02, 1.0, 0.102),
            loc="lower left",
            ncol=3,
            mode="expand",
            borderaxespad=0.0,
            fontsize=fs,
        )

    # FDC axes[2]
    for label, color in zip(labels, colors):
        axes[2].plot(
            np.arange(0, len(dsq.time)) / (len(dsq.time) + 1),
            dsq.Q.sel(runs=label).sortby(
                dsq.Q.sel(
                    runs=label,
                ),
                ascending=False,
            ),
            color=color,
            linestyle="-",
            linewidth=lw,
            label=label,
        )

    if not skip_obs:
        axes[2].plot(
            np.arange(0, len(dsq.time)) / (len(dsq.time) + 1),
            dsq.Q.sel(runs="Obs.").sortby(
                dsq.Q.sel(
                    runs="Obs.",
                ),
                ascending=False,
            ),
            color="k",
            linestyle=":",
            linewidth=lw,
            label="Obs.",
        )
    axes[2].set_xlabel("Exceedence probability (-)", fontsize=fs)
    axes[2].set_ylabel("Q (m$^3$s$^{-1}$)", fontsize=fs)

    # FDClog axes[3]
    for label, color in zip(labels, colors):
        axes[3].plot(
            np.arange(0, len(dsq.time)) / (len(dsq.time) + 1),
            np.log(
                dsq.Q.sel(runs=label).sortby(
                    dsq.Q.sel(
                        runs=label,
                    ),
                    ascending=False,
                )
            ),
            color=color,
            linestyle="-",
            linewidth=lw,
            label=label,
        )
    if not skip_obs:
        axes[3].plot(
            np.arange(0, len(dsq.time)) / (len(dsq.time) + 1),
            np.log(
                dsq.Q.sel(runs="Obs.").sortby(
                    dsq.Q.sel(
                        runs="Obs.",
                    ),
                    ascending=False,
                )
            ),
            color="k",
            linestyle=":",
            linewidth=lw,
            label="Obs.",
        )
    axes[3].set_xlabel("Exceedence probability (-)", fontsize=fs)
    axes[3].set_ylabel("log(Q)", fontsize=fs)

    # max annual axes[4]
    time_slice = slice(
        f"{str(dsq['time.year'][0].values)}-09-01",
        f"{str(dsq['time.year'][-1].values)}-08-31",
    )
    if (dsq['time.year'][-1].values - dsq['time.year'][-1].values) == 0:
        stime = np.datetime_as_string(dsq["time"][0].values, unit="D")
        etime = np.datetime_as_string(dsq["time"][-1].values, unit="D")
        time_slice = slice(
            stime,
            etime,
        )
    dsq_max = (
        dsq.sel(
            time=time_slice,
        )
        .resample(time="AS-Sep")
        .max("time")
    )
    if not skip_obs:
        for label, color in zip(labels, colors):
            axes[4].plot(
                dsq_max.Q.sel(runs="Obs."),
                dsq_max.Q.sel(runs=label),
                color=color,
                marker="o",
                linestyle="None",
                linewidth=lw,
                label=label,
            )
    axes[4].plot(
        [0, max_y * 1.1], [0, max_y * 1.1], color="0.5", linestyle="--", linewidth=1
    )
    axes[4].set_xlim([0, max_y * 1.1])
    axes[4].set_ylim([0, max_y * 1.1])
    # R2 score add!
    if len(labels) + 1 < legend_threshold and not skip_obs:
        text = ""
        for label in labels:
            r2_score = rsquared(
                dsq_max["Q"].sel(runs="Obs."), dsq_max["Q"].sel(runs=label)
            )
            text += f"R$_2$ {label} = {r2_score:.2f} \n"
        axes[4].text(
            0.01,
            0.99,
            text,
            fontsize=fs,
            ha="left",
            va="top",
            transform=axes[4].transAxes,
        )

    # # add MHQ
    # mhq = dsq_max.mean("time")
    # if not skip_obs:
    #     for label, color in zip(labels, colors):
    #         axes[4].plot(
    #             mhq.Q.sel(runs="Obs."),
    #             mhq.Q.sel(runs=label),
    #             color="black",
    #             marker=">",
    #             linestyle="None",
    #             linewidth=lw,
    #             label=label,
    #             markersize=6,
    #         )
    # axes[4].plot(mhq.Q.sel(runs = 'Obs.'), mhq.Q.sel(runs = label_01), color = 'grey', marker = '^', linestyle = 'None', linewidth = lw, label = label_01, markersize = 6)
    # labels
    axes[4].set_ylabel("Sim. max annual Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[4].set_xlabel("Obs. max annual Q (m$^3$s$^{-1}$)", fontsize=fs)

    # nm7q axes[5]
    dsq_nm7q = dsq.rolling(time=window).mean().resample(time="A").min("time")
    max_ylow = dsq_nm7q["Q"].max().values
    if not skip_obs:
        for label, color in zip(labels, colors):
            axes[5].plot(
                dsq_nm7q.Q.sel(runs="Obs."),
                dsq_nm7q.Q.sel(runs=label),
                color=color,
                marker="o",
                linestyle="None",
                linewidth=lw,
                label=label,
            )
    # axes[5].plot(dsq_nm7q.Q.sel(runs = 'Obs.'), dsq_nm7q.Q.sel(runs = label_01), color = color_01, marker = '.', linestyle = 'None', linewidth = lw, label = label_01)
    axes[5].plot(
        [0, max_ylow * 1.1],
        [0, max_ylow * 1.1],
        color="0.5",
        linestyle="--",
        linewidth=1,
    )
    axes[5].set_xlim([0, max_ylow * 1.1])
    axes[5].set_ylim([0, max_ylow * 1.1])
    # R2 score add!
    if len(labels) + 1 < legend_threshold and not skip_obs:
        text = ""
        for label in labels:
            r2_score = rsquared(
                dsq_nm7q["Q"].sel(runs="Obs."), dsq_nm7q["Q"].sel(runs=label)
            )
            text += f"R$_2$ {label} = {r2_score:.2f} \n"
        axes[5].text(
            0.01,
            0.99,
            text,
            fontsize=fs,
            ha="left",
            va="top",
            transform=axes[5].transAxes,
        )
    # labels
    axes[5].set_ylabel("Simulated NM7Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[5].set_xlabel("Observed NM7Q (m$^3$s$^{-1}$)", fontsize=fs)

    # gumbel high axes[6]
    a = 0.3
    b = 1.0 - 2.0 * a
    ymin, ymax = 0, max_y
    p1 = ((np.arange(1, len(dsq_max.time) + 1.0) - a)) / (len(dsq_max.time) + b)
    RP1 = 1 / (1 - p1)
    gumbel_p1 = -np.log(-np.log(1.0 - 1.0 / RP1))
    ts = [2.0, 5.0, 10.0, 30.0]  # ,30.,100.,300.,1000.,3000.,10000.,30000.]
    # plot
    for label, color in zip(labels, colors):
        axes[6].plot(
            gumbel_p1,
            dsq_max["Q"].sel(runs=label).sortby(dsq_max["Q"].sel(runs=label)),
            marker="o",
            color=color,
            linestyle="None",
            label=label,
            markersize=4,
        )
    if not skip_obs:
        axes[6].plot(
            gumbel_p1,
            dsq_max["Q"].sel(runs="Obs.").sortby(dsq_max["Q"].sel(runs="Obs.")),
            marker="+",
            color="k",
            linestyle="None",
            label="Obs.",
            markersize=6,
        )
    # axes[6].plot(gumbel_p1, dsq_max['Q'].sel(runs = label_01).sortby(dsq_max['Q'].sel(runs = label_01)), marker = '.', color = color_01, linestyle = 'None', label = label_01, markersize = 3)

    for t in ts:
        axes[6].axvline(-np.log(-np.log(1 - 1.0 / t)), c="0.5", alpha=0.4)
        axes[6].text(
            -np.log(-np.log(1 - 1.0 / t)),
            ymax * 1.05,
            f"{t:.0f}y",
            fontsize=fs,
            ha="center",
            va="bottom",
        )

    axes[6].set_ylabel("max. annual Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[6].set_xlabel("Plotting position and associated return period", fontsize=fs)

    # gumbel low axes[7]
    a = 0.3
    b = 1.0 - 2.0 * a
    _, ymax = 0, max_ylow
    p1 = ((np.arange(1, len(dsq_nm7q.time) + 1.0) - a)) / (len(dsq_nm7q.time) + b)
    RP1 = 1 / (1 - p1)
    gumbel_p1 = -np.log(-np.log(1.0 - 1.0 / RP1))
    ts = [2.0, 5.0, 10.0, 30.0]  # ,30.,100.,300.,1000.,3000.,10000.,30000.]
    # plot
    for label, color in zip(labels, colors):
        axes[7].plot(
            gumbel_p1,
            dsq_nm7q["Q"]
            .sel(runs=label)
            .sortby(dsq_nm7q["Q"].sel(runs=label), ascending=False),
            marker="o",
            color=color,
            linestyle="None",
            label=label,
            markersize=4,
        )
    if not skip_obs:
        axes[7].plot(
            gumbel_p1,
            dsq_nm7q["Q"]
            .sel(runs="Obs.")
            .sortby(dsq_nm7q["Q"].sel(runs="Obs."), ascending=False),
            marker="+",
            color="k",
            linestyle="None",
            label="Obs.",
            markersize=6,
        )
    # axes[7].plot(gumbel_p1, dsq_nm7q['Q'].sel(runs = label_01).sortby(dsq_nm7q['Q'].sel(runs = label_01), ascending=False), marker = '.', color = color_01, linestyle = 'None', label = label_01, markersize = 3)

    for t in ts:
        axes[7].axvline(-np.log(-np.log(1 - 1.0 / t)), c="0.5", alpha=0.4)
        axes[7].text(
            -np.log(-np.log(1 - 1.0 / t)),
            ymax * 1.05,
            f"{t:.0f}y",
            fontsize=fs,
            ha="center",
            va="bottom",
        )

    axes[7].set_ylabel("NM7Q (m$^3$s$^{-1}$)", fontsize=fs)
    axes[7].set_xlabel("Plotting position and associated return period", fontsize=fs)

    # cum axes[8]
    if not skip_obs:
        dsq["Q"].sel(runs="Obs.").cumsum("time").plot(
            ax=axes[8], color="k", linestyle=":", linewidth=lw, label="Obs."
        )
    for label, color in zip(labels, colors):
        dsq["Q"].sel(runs=label).cumsum("time").plot(
            ax=axes[8], color=color, linestyle="-", linewidth=lw, label=label
        )
    # dsq['Q'].sel(runs = label_01).cumsum('time').plot(ax=axes[8], color = color_01, linestyle = '--', linewidth = lw, label = label_01)
    axes[8].set_xlabel("")
    axes[8].set_ylabel("Cum. Q (m$^3$s$^{-1}$)", fontsize=fs)

    # performance measures KGE, NSE, NSlogQ, NSlog_mm7q axes[9]
    #     sns.boxplot(ax=axes[9], data = df_perf, x = 'metrics', hue = 'runs', y = 'performance')
    offsets = np.linspace(-0.25, 0.25, len(labels))
    # kge
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            0.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="KGE")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    # axes[9].plot(5.2, dsq['performance'].loc[dict(runs = label_01, metrics = 'KGE')], color = color_01, marker = 'o', linestyle = 'None', linewidth = lw, label = label_01)
    # nse
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            2.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="NSE")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    # axes[9].plot(1.2, dsq['performance'].loc[dict(runs = label_01, metrics = 'NSE')], color = color_01, marker = 'o', linestyle = 'None', linewidth = lw, label = label_01)
    # nselog
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            4.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="NSElog")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    # axes[9].plot(3.2, dsq['performance'].loc[dict(runs = label_01, metrics = 'NSElog')], color = color_01, marker = 'o', linestyle = 'None', linewidth = lw, label = label_01)
    # nselog_mm7q
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            6.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="NSElog_mm7q")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    axes[9].set_xticks([1, 3, 5, 7])
    axes[9].set_xticklabels(["KGE", "NSE", "NSElog", "NSElog_mm7q"])
    axes[9].set_ylim([0, 1])
    axes[9].set_ylabel("Performance", fontsize=fs)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=fs)
        ax.set_title("")
        ax.yaxis.offsetText.set_fontsize(fs)

    # plt.tight_layout()
    if save:
        fig.savefig(
            os.path.join(Folder_out, f"signatures_{station_name}_{station_id}.png"),
            dpi=300,
        )
        # fig.close()
    return metrics_df
def save_metrics_to_excel(all_metrics_df, output_dir):
    excel_path = os.path.join(output_dir, 'all_station_metrics.xlsx')
    all_metrics_df.to_excel(excel_path, index=False)
    print(f"All metrics saved to {excel_path}")

def plot_hydro(
    dsq,
    start_long,
    end_long,
    start_1,
    end_1,
    start_2,
    end_2,
    start_3,
    end_3,
    labels,
    colors,
    Folder_out,
    station_name,
    station_id,
    lw=0.8,
    fs=8,
    save=False,
):
    # fig, axes = plt.subplots(4, 1, figsize=(16 / 2.54, 20 / 2.54))
    dsq, labels, colors = process_labels_and_colors(dsq, labels, colors)
    
    fig = plt.figure(
        "timeseries", clear=True, tight_layout=True, figsize=(16 / 2.54, 20 / 2.54)
    )
    axes = fig.subplots(4, 1)

    fig.suptitle(station_name)

    # period long
    dsq["Q"].sel(runs="Obs.", time=slice(start_long, end_long)).plot(
        ax=axes[0], label="Obs.", linewidth=lw, color="k", linestyle="--"
    )
    # period zoom 1
    dsq["Q"].sel(runs="Obs.", time=slice(start_1, end_1)).plot(
        ax=axes[1], label="Obs.", linewidth=lw, color="k", linestyle="--"
    )
    # period zoom 2
    dsq["Q"].sel(runs="Obs.", time=slice(start_2, end_2)).plot(
        ax=axes[2], label="Obs.", linewidth=lw, color="k", linestyle="--"
    )
    # period zoom 3
    dsq["Q"].sel(runs="Obs.", time=slice(start_3, end_3)).plot(
        ax=axes[3], label="Obs.", linewidth=lw, color="k", linestyle="--"
    )

    # for label, color in zip(labels, colors):
    for i, label in enumerate(labels):
        color = colors[i]
        dsq["Q"].sel(runs=label, time=slice(start_long, end_long)).plot(
            ax=axes[0], label=label, linewidth=lw, color=color
        )
        dsq["Q"].sel(runs=label, time=slice(start_1, end_1)).plot(
            ax=axes[1], label=label, linewidth=lw, color=color
        )
        dsq["Q"].sel(runs=label, time=slice(start_2, end_2)).plot(
            ax=axes[2], label=label, linewidth=lw, color=color
        )
        dsq["Q"].sel(runs=label, time=slice(start_3, end_3)).plot(
            ax=axes[3], label=label, linewidth=lw, color=color
        )

    axes[0].set_xlim(
        datetime.strptime(start_long, "%Y-%m-%d"),
        datetime.strptime(end_long, "%Y-%m-%d"),
    )
    axes[1].set_xlim(
        datetime.strptime(start_1, "%Y-%m-%d"), datetime.strptime(end_1, "%Y-%m-%d")
    )
    axes[2].set_xlim(
        datetime.strptime(start_2, "%Y-%m-%d"), datetime.strptime(end_2, "%Y-%m-%d")
    )
    axes[3].set_xlim(
        datetime.strptime(start_3, "%Y-%m-%d"), datetime.strptime(end_3, "%Y-%m-%d")
    )

    for ax in axes:
        ax.tick_params(axis="both", labelsize=fs)
        ax.set_ylabel("Q (m$^3$s$^{-1}$)", fontsize=fs)
        ax.set_xlabel("", fontsize=fs)
        ax.set_title("")

    axes[0].legend(fontsize=fs)
    # plt.tight_layout()

    if save:
        station_name = station_name.replace(" ", "_")
        fig.savefig(os.path.join(Folder_out, f"hydro_{station_name}_{station_id}.png"), dpi=300)
        # fig.close()


if __name__ == "__main__":

    import platform
    if platform.system() == "Windows":
        DRIVE = "p:/"
        PLATFORM = "Windows"
    elif platform.system() == "Linux":
        DRIVE = "/p"
        PLATFORM = "Linux"
        
    work_dir = Path(DRIVE, '11209265-grade2023', 'wflow', 'RWSOS_Calibration', 'meuse_random_spider')
    ds_fn = work_dir / 'data/4-output/combined_output.nc'
    output_dir = work_dir / 'data/5-visualization/signature'
    output_dir.mkdir(parents=True, exist_ok=True)
    GaugeToPlot = Path(work_dir,'data', '4-output','wflow_id_add_HBV_new.csv')
    
    ds=xr.open_dataset(ds_fn)
    gauges = pd.read_csv(GaugeToPlot)['wflow_id'].tolist()
    
    # Colorblind-friendly colors
    colorblind_friendly_colors = [
        "#ff7f00",  # HBV -- RWSOS report colot
        "#b2df8a",  # Base -- RWSOS report color
        "#56B4E9",  # Sky Blue vv colorblind pallette vv
        "#009E73",  # Bluish Green
        "#F0E442",  # Yellow
        "#0072B2",  # Blue
        "#D55E00",  # Vermilion
        "#CC79A7",  # Reddish Purple
        "#999999",  # Gray
        "#44AA99",  # Teal
        "#661100",  # Brown
        "#6699CC",  # Light Blue ^^ colorblind pallette ^^
        "#1f78b4",  # Front facing -- RWSOS report color (green)
    ]
    labels = ["HBV", "base_model", *[f"Top_{i+1}" for i in range(10)]]
    lc = {label: color for label, color in zip(labels, colorblind_friendly_colors)}
    
    # Use this line to replace the original colors with the colorblind-friendly ones
    colors = colorblind_friendly_colors
    plot_colors = colors[:]
    start_year = 2006
    end_year = 2016
    
    start_long = f"{start_year}-01-01"
    end_long = f"{end_year}-01-01"
    start_1 = f"{start_year}-08-01"
    end_1 = f"{start_year+1}-07-30"
    start_2 = f"{2011}-08-01"
    end_2 = f"{2012}-07-30"
    start_3 = f"{2014}-08-01"
    end_3 = f"{2015}-07-30"
    
    ds = ds.sel(time=slice(f"{start_year}-01-01", f"{end_year}-01-01"))
    start_time = time.time()
    all_metrics_df = pd.DataFrame()
    for station_id in tqdm(gauges, desc="Processing gauges", total=len(gauges), colour='green'):
        # Get data for that station
        dsq = ds.sel(wflow_id=station_id)
        station_name = pd.read_csv(GaugeToPlot).set_index('wflow_id')['location'].loc[station_id]
        # Prep labels
        labels = {"HBV": "HBV", "base_model": "2023 Base", "Top_10": "Calibrated"}
        plot_colors = {label: lc[label] for label in labels.keys()}
        
        metrics_df = plot_signatures(
            dsq=dsq,
            labels=labels,
            colors=plot_colors,
            Folder_out=output_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
        )
        all_metrics_df = pd.concat([all_metrics_df, metrics_df], ignore_index=True)
        
        plot_hydro(
            dsq=dsq,
            labels=labels,
            colors=plot_colors,
            Folder_out=output_dir,
            station_name=station_name,
            station_id=station_id,
            save=True,
            start_long=start_long,
            end_long=end_long,
            start_1=start_1,
            end_1=end_1,
            start_2=start_2,
            end_2=end_2,
            start_3=start_3,
            end_3=end_3,
        )
        print(f"Finished processing {station_name}")
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2f} seconds")


