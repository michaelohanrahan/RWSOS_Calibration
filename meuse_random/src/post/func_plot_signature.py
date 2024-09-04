import os
from datetime import datetime

import hydromt
import numpy as np
import matplotlib.pyplot as plt
import scipy
import xarray as xr


def rsquared(x, y):
    """Return R^2 where x and y are array-like."""

    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(x, y)
    return r_value**2


def plot_signatures(
    dsq,
    labels,
    colors,
    Folder_out,
    station_name,
    station_id,
    lw=0.8,
    fs=8,
    save=False,
    legend_threshold=5,
):
    if xr.infer_freq(dsq.time).lower() == "D":
        window = 7
    elif xr.infer_freq(dsq.time).lower() == "h":
        window = 7 * 24
    elif xr.infer_freq(dsq.time).lower() == "3h":
        window = 7 * int(24 / 3)

    # first calc some signatures
    dsq["metrics"] = ["KGE", "NSE", "NSElog", "RMSE", "MSE"]
    dsq["performance"] = (
        ("runs", "metrics"),
        np.zeros((len(dsq.runs), len(dsq.metrics))) * np.nan,
    )

    # dsplot = dsq.copy(deep=True)
    _tmp = dsq.sel(time=dsq.time[~dsq.Q.sel(runs="Obs.").isnull()].values)
    if len(_tmp.time) == 0:
        skip_obs = True

    else:
        skip_obs = False
        dsq = dsq.sel(time=dsq.time[~dsq.Q.sel(runs="Obs.").isnull()].values)

    # perf metrics for single station
    if not skip_obs:
        for label in labels:
            # nse
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
        #         print(nse.values, nselog.values, kge['kge'].values, rmse.values, mse.values)

    # needed later for sns boxplot
    #     df_perf = pd.DataFrame()
    #     for label in [label_00, label_01]:
    #         df = dsq['performance'].sel(runs = label, metrics = ['NSE', 'NSElog', 'KGE']).to_dataframe()
    #         df_perf = pd.concat([df,df_perf])

    fig = plt.figure(
        "signatures", clear=True, figsize=(16 / 2.54, 22 / 2.54), tight_layout=True
    )
    # fig, axes = plt.subplots(5, 2, figsize=(16 / 2.54, 22 / 2.54))
    axes = fig.subplots(5, 2)
    axes = axes.flatten()

    fig.suptitle(station_name)

    # daily against each other axes[0]
    if not skip_obs:
        for label, color in zip(labels, colors):
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

    # performance measures NS, NSlogQ, KGE, axes[9]
    #     sns.boxplot(ax=axes[9], data = df_perf, x = 'metrics', hue = 'runs', y = 'performance')
    # nse
    offsets = np.linspace(-0.25, 0.25, len(labels))
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            0.8 + offset,
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
            2.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="NSElog")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    # axes[9].plot(3.2, dsq['performance'].loc[dict(runs = label_01, metrics = 'NSElog')], color = color_01, marker = 'o', linestyle = 'None', linewidth = lw, label = label_01)
    # kge
    for label, color, offset in zip(labels, colors, offsets):
        axes[9].plot(
            4.8 + offset,
            dsq["performance"].loc[dict(runs=label, metrics="KGE")],
            color=color,
            marker="o",
            linestyle="None",
            linewidth=lw,
            label=label,
        )
    # axes[9].plot(5.2, dsq['performance'].loc[dict(runs = label_01, metrics = 'KGE')], color = color_01, marker = 'o', linestyle = 'None', linewidth = lw, label = label_01)
    axes[9].set_xticks([1, 3, 5])
    axes[9].set_xticklabels(["NSE", "NSElog", "KGE"])
    axes[9].set_ylim([0, 1])
    axes[9].set_ylabel("Performance", fontsize=fs)

    for ax in axes:
        ax.tick_params(axis="both", labelsize=fs)
        ax.set_title("")
        ax.yaxis.offsetText.set_fontsize(fs)

    # plt.tight_layout()
    if save:
        fig.savefig(
            os.path.join(Folder_out, f"signatures_{station_id}.png"),
            dpi=300,
        )
        # fig.close()


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
        fig.savefig(os.path.join(Folder_out, f"hydro_{station_id}.png"), dpi=300)
        # fig.close()
