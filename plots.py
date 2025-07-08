#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 14:45:47 2021

@author: ja17375

Functions to plot shear-wave ata and shear-wave splitting results
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import obspy


def diagnostic_plot(st_in, st_corr_in, result, event_time):
    """
    Produces a diagnostic plot for each shear-wave splitting measurment

    Shows:
        - corrected/uncorrected traces, with windows and S pick
        - corrected/uncorrected particle motions
        - normalised egigenvalue surface

    Parameters:
        st : Stream (Obspy)
            input data to SHEBA
        st_corr : Stream (Obspy)
            output data corrected by SHEBA for splitting
        result : Dataset (netCDF4)
            dataset netCDF4 files pre read in
        event_time :
    """

    plt.close()

    # trim st, st_corr to 10% of window before and after
    st = st_in.copy()
    st.trim(event_time + result.wbeg * 0.9, event_time + result.wend * 1.1)
    st_corr = st_corr_in.copy()
    st_corr.trim(event_time + result.wbeg * 0.9, event_time + result.wend * 1.1)

    fig = plt.figure(layout="constrained", figsize=(13, 9))
    gs = GridSpec(6, 6, figure=fig)
    # Input data ZNE
    ax1 = fig.add_subplot(gs[0, 0:3])
    _plot_traces(
        st,
        show_final_window=True,
        axes=ax1,
        event_time=event_time,
        cmp_orientation="NE",
        wbeg=result.wbeg,
        wend=result.wend,
    )
    ax1.set_xlim([result.wbeg * 0.95, result.wend * 1.05])
    time_str = event_time.strftime("%Y-%m-%d %H:%M:%S %Z")
    ax1.set_title(
        f"Input S. Event: {time_str}" + f" Station: {result.station.strip()}",
        fontsize=12,
    )
    # Input data RT
    ax2 = fig.add_subplot(gs[1, 0:3], sharey=ax1, sharex=ax1)
    _plot_traces(
        st,
        show_final_window=True,
        axes=ax2,
        event_time=event_time,
        cmp_orientation="RT",
        spol=result.spol,
        wbeg=result.wbeg,
        wend=result.wend,
    )
    ax2.set_title(f"Input S. Radial-Transverse. Spol: {result.spol:4.2f}°")
    # Corrected data ZNE
    ax3 = fig.add_subplot(gs[0, 3:], sharey=ax1, sharex=ax1)
    _plot_traces(
        st_corr,
        show_final_window=True,
        axes=ax3,
        event_time=event_time,
        cmp_orientation="NE",
        wbeg=result.wbeg,
        wend=result.wend,
    )
    fast_res = rf"$\phi_f = {result.fast:4.2f}\pm{result.dfast:4.2f}$°"
    dt_res = rf"$\delta t = {result.tlag:4.3f}\pm{result.dtlag:4.3f}$ s"
    ax3.set_title(f"Corrected S. {fast_res}, {dt_res}")
    # Corrected data RT
    ax4 = fig.add_subplot(gs[1, 3:], sharey=ax1, sharex=ax1)
    _plot_traces(
        st_corr,
        show_final_window=True,
        axes=ax4,
        event_time=event_time,
        cmp_orientation="RT",
        spol=result.spol,
        wbeg=result.wbeg,
        wend=result.wend,
    )

    ax4.set_title(f"Corrected S. Radial-Transverse. Spol: {result.spol:4.2f}°")

    # Fast-slow orig
    ax5 = fig.add_subplot(gs[2, 0])
    _plot_fast_slow(
        st,
        event_time,
        phi=result.fast,
        axes=ax5,
        norm=True,
        wbeg=result.wbeg,
        wend=result.wend,
    )
    # Fast-slow corrected normalised
    ax6 = fig.add_subplot(gs[2, 1])
    _plot_fast_slow(
        st_corr,
        event_time,
        phi=result.fast,
        tshift=result.tlag,
        norm=True,
        axes=ax6,
        wbeg=result.wbeg,
        wend=result.wend,
    )
    # Fast Slow corrected not normalised
    ax7 = fig.add_subplot(gs[2, 2])
    _plot_fast_slow(
        st_corr,
        event_time,
        phi=result.fast,
        axes=ax7,
        wbeg=result.wbeg,
        wend=result.wend,
    )
    # Particle motion uncorrected
    ax8 = fig.add_subplot(gs[3, 0])
    _ppm(ax8, st, event_time, wbeg=result.wbeg, wend=result.wend)

    # Particle motion corrected
    ax9 = fig.add_subplot(gs[3, 1], sharex=ax8, sharey=ax8)
    _ppm(ax9, st_corr, event_time, wbeg=result.wbeg, wend=result.wend)

    # Eigenvalue grid surface

    ax10 = fig.add_subplot(gs[2:, 3:6])
    phis = result.variables["fast_vector"][:]
    dts = result.variables["tlag_vector"][:]
    PHI, TLAG = np.meshgrid(phis, dts)

    C1 = ax10.contour(
        TLAG,
        PHI,
        result.variables["lam2_norm_grid"][:],
        levels=[2, 4, 6, 8, 10, 12, 16, 20],
        colors="black",
    )
    ax10.clabel(C1, fontsize=8, fmt="%2.0f")
    C2 = ax10.contour(
        TLAG, PHI, result.variables["lam2_norm_grid"][:], levels=[1], linewidths=[3]
    )
    ax10.clabel(C2, fontsize=10, fmt="%2.0f")
    ax10.plot(result.tlag, result.fast, "x", color="royalblue")
    ax10.set_ylabel("Fast polarisation [°]")
    ax10.set_xlabel("Delay time [s]")
    ax10.set_xlim([dts.min(), dts.max()])
    ax10.set_ylim([phis.min(), phis.max()])

    # Window clustser fast
    # shift to start from 1
    wind_num = np.arange(0, result.dimensions["window"].size) + 1
    ax11 = fig.add_subplot(gs[4, 0:3])
    ax11.errorbar(
        x=wind_num,
        y=result.variables["mw_fast"][:],
        yerr=result.variables["mw_dfast"][:],
    )
    ax11.errorbar(
        x=wind_num[result.best_window - 1],
        y=result.variables["mw_fast"][result.best_window - 1],
        yerr=result.variables["mw_dfast"][result.best_window - 1],
        color="tab:red",
    )
    ax11.set_ylim([phis.min(), phis.max()])
    ax11.set_ylabel(r"$\phi_f$ [°]")
    # dt for all windows
    ax12 = fig.add_subplot(gs[5, 0:3])
    ax12.errorbar(
        x=wind_num,
        y=result.variables["mw_tlag"][:],
        yerr=result.variables["mw_dtlag"][:],
    )
    ax12.errorbar(
        wind_num[result.best_window - 1],
        y=result.variables["mw_tlag"][result.best_window - 1],
        yerr=result.variables["mw_dtlag"][result.best_window - 1],
        color="tab:red",
    )
    ax12.set_ylim([dts.min(), dts.max()])
    ax12.set_ylabel(r"$\delta t$ [s]")
    ax12.set_xlabel("Window #")
    # Clusters
    ax13 = fig.add_subplot(gs[3, 2])
    ax13.scatter(
        x=result.variables["cluster_xc0"][:], y=result.variables["cluster_yc0"][:]
    )
    ax13.plot(
        result.variables["cluster_xc0"][result.best_cluster - 1],
        result.variables["cluster_yc0"][result.best_cluster - 1],
        "x",
        color="red",
    )
    ax13.set_xlim([dts.min(), dts.max()])
    ax13.set_ylim([phis.min(), phis.max()])
    ax13.set_ylabel(r"$\delta t$ [s]")
    ax13.set_ylabel(r"$\phi_f$ [°]")

    return fig


def _plot_fast_slow(st, event_time, phi, norm=False, **kwargs):
    """
    function to plot fast and slow shear waves in
    window
    """
    if "axes" not in kwargs:
        # no axes provided so make our own
        fig, ax = plt.subplots(1, 1)
    else:
        ax = kwargs["axes"]
    st_plot = st.copy()
    if "wbeg" in kwargs:
        st_plot.trim(event_time + kwargs["wbeg"], event_time + kwargs["wend"])
        ax.set_xlim([kwargs["wbeg"], kwargs["wend"]])
    else:
        st_plot.trim(
            event_time + st[0].stats.sac["a"], event_time + st[0].stats.sac["f"]
        )
    if norm:
        st_plot.normalize()
    # rotate to fast-slow
    st_plot.rotate(method="NE->RT", back_azimuth=(phi + 360) % 360)

    times = st_plot[0].times(reftime=event_time)
    ax.plot(times, st_plot[0].data, label="S1", color="black", linestyle="-")
    ax.plot(times, st_plot[1].data, label="S2", color="black", linestyle="--")
    # set axis label
    if "units" not in kwargs:
        kwargs["units"] = "s"
    ax.set_xlabel(f'Time ({kwargs["units"]})')

    return


def _plot_traces(st, event_time, cmp_orientation="NE", **kwargs):
    """
    function to plot shear-wave traces

    Parameters:
    ----------
    st :
        obspy Stream conatining waveform data to plot
    """
    if "axes" not in kwargs:
        # no axes provided so make our own
        fig, ax = plt.subplots(1, 1)
    else:
        ax = kwargs["axes"]

    if cmp_orientation == "RT":
        st_plot = st.copy()
        if "spol" in kwargs:
            st_plot.rotate(method="NE->RT", back_azimuth=kwargs["spol"])
        else:
            st_plot.rotate(method="NE->RT", back_azimuth=st[0].stats.sac["baz"])
    else:
        st_plot = st.copy()

    st_plot.trim(
        event_time + 0.95 * kwargs["wbeg"],
        event_time + 1.05 * kwargs["wend"],
    )

    times = st_plot[0].times(reftime=event_time)
    ax.plot(times, st_plot[0].data, label=st_plot[0].stats.channel, color="dodgerblue")
    ax.plot(times, st_plot[1].data, label=st_plot[1].stats.channel, color="darkorange")
    # trim to window

    # set y axis limits
    ylim = 1.1 * np.max(
        [
            st_plot[0].data.max(),
            st_plot[1].data.max(),
            abs(st_plot[0].data.min()),
            abs(st_plot[1].data.min()),
        ]
    )
    ax.set_ylim([-ylim, ylim])

    if "show_window_range" in kwargs:
        # window range should be a list of [wbeg1 wend1 wbeg2 wend2]
        for marker in ["user0", "user1", "user2", "user3"]:
            ax.axvline(
                x=st_plot[0].stats.sac[marker],
                linewidth=1,
                color="black",
                linestyle="--",
            )

    if "show_final_window" in kwargs:
        try:
            ax.axvline(
                x=st_plot[0].stats.sac["a"], linewidth=1, color="black", linestyle="-"
            )
            ax.axvline(
                x=st_plot[0].stats.sac["f"], linewidth=1, color="black", linestyle="-"
            )
        except KeyError:
            print("No final window in SAC file")

    # set axis label
    if "units" not in kwargs:
        kwargs["units"] = "s"

    ax.set_xlabel(f'Time relative to origin ({kwargs["units"]})')
    ax.legend(framealpha=0.75, loc="upper right", ncols=2)

    return


def _ppm(ax, st, event_time, **kwargs):
    st_plot = st.copy()
    st_plot.normalize()

    if "wbeg" in kwargs:
        st_plot.trim(event_time + kwargs["wbeg"], event_time + kwargs["wend"])
    else:
        st_plot.trim(
            event_time + st[0].stats.sac["a"], event_time + st[0].stats.sac["f"]
        )
    trN = st_plot.select(channel="??N")[0]
    trE = st_plot.select(channel="??E")[0]
    ax.plot(trE.data, trN.data)
    ax.set_xlabel("East")
    ax.set_ylabel("North")
    # make x and y axes the same extent
    axis_extent = 1.05 * max(
        abs(trE.data.min()), abs(trN.data.min()), trE.data.max(), trN.data.max()
    )
    ax.set_xlim([-axis_extent, axis_extent])
    ax.set_ylim([-axis_extent, axis_extent])

    return


if __name__ == "__main__":
    st = obspy.read(
        "/Users/eart0593/Projects/SHARP/splitting/data/ML_gte_0/*/GB_MONM_20180420145936.*"
    )

    _plot_traces(st)
    plt.show()
