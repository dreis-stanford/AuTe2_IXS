"""2D visualization of phonon frequencies and IXS intensity across the
BL43LXU analyzer array.

For each phonon branch, plots a contour map of the mode frequency over the
analyzer array's (delta_tth, delta_gamma) layout -- the angular offsets of
each analyzer from the array center, computed by
SixCircleInterface.analyzer_array_offsets(). The IXS Stokes cross-section at
each analyzer is overlaid as an open circle whose area encodes
log10(IXS_stokes), so both quantities can be read off the same plot.
"""
import numpy as np


def plot_analyzer_array(array_results, freq_unit='meV', Q_label=None, grid=None):
    """
    Build a grid figure with one frequency-contour subplot per phonon branch.

    Parameters
    ----------
    array_results : list of dict
        Output of SingleQAnalyzer.analyze_array(). Each entry must have
        'dtth', 'dgam', 'xi', 'zi' (from analyzer_array_offsets), the
        frequency array for `freq_unit`, and 'IXS_stokes'.
    freq_unit : str
        'meV', 'THz', or 'cm-1' -- selects frequencies_<unit> and axis labels.
    Q_label : str, optional
        Text describing the center Q, shown in the figure title.
    grid : dict, optional
        Output of SingleQAnalyzer.analyze_array_grid() -- a denser
        (dtth, dgam) grid spanning a wider region than the analyzer array
        itself. If given, the frequency contours are drawn over this wider
        grid instead of just the 28 analyzer positions, with the actual
        analyzers overlaid as circles as usual.

    Returns
    -------
    fig : matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    freq_key = {'meV': 'frequencies_meV',
                 'cm-1': 'frequencies_cm',
                 'THz': 'frequencies_THz'}[freq_unit]
    freq_label = {'meV': 'meV', 'cm-1': r'cm$^{-1}$', 'THz': 'THz'}[freq_unit]

    n_modes = len(array_results[0][freq_key])
    nx = max(r['xi'] for r in array_results) + 1
    nz = max(r['zi'] for r in array_results) + 1

    dtth = np.full((nz, nx), np.nan)
    dgam = np.full((nz, nx), np.nan)
    freq = np.full((nz, nx, n_modes), np.nan)
    ixs = np.full((nz, nx, n_modes), np.nan)

    for r in array_results:
        xi, zi = r['xi'], r['zi']
        dtth[zi, xi] = r['dtth']
        dgam[zi, xi] = r['dgam']
        freq[zi, xi, :] = r[freq_key]
        ixs[zi, xi, :] = r['IXS_stokes']

    # log10(IXS) -> marker area, shared scale across all modes so sizes are
    # comparable subplot-to-subplot. Some modes are selection-rule-forbidden
    # at this Q (IXS ~ 1e-30, numerically zero), so the true minimum is not a
    # useful floor -- it would compress all the "real" signal (which spans
    # ~4-5 decades) into a sliver near the top of the size range. Instead,
    # cap the dynamic range to a fixed number of decades below the brightest
    # point; anything dimmer than that floors out at the minimum marker size.
    DYNAMIC_RANGE_DECADES = 4.0
    log_ixs = np.log10(np.clip(ixs, 1e-300, None))
    finite = log_ixs[np.isfinite(log_ixs)]
    vmax = finite.max() if finite.size else 0.0
    vmin = vmax - DYNAMIC_RANGE_DECADES

    def sizes_for(log_vals):
        frac = np.clip((log_vals - vmin) / (vmax - vmin), 0.0, 1.0)
        return 15.0 + frac * 200.0  # marker area, points^2

    ncols = 3
    nrows = int(np.ceil(n_modes / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.0 * ncols, 2.6 * nrows),
                              squeeze=False)

    if grid is not None:
        grid_dtth, grid_dgam, grid_freq = grid['dtth'], grid['dgam'], grid[freq_key]

    for imode in range(n_modes):
        ax = axes[imode // ncols, imode % ncols]
        if grid is not None:
            cf = ax.contourf(grid_dtth, grid_dgam, grid_freq[:, :, imode],
                              levels=12, cmap='viridis')
        else:
            cf = ax.contourf(dtth, dgam, freq[:, :, imode], levels=12, cmap='viridis')
        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label(freq_label)

        ax.scatter(dtth, dgam, s=sizes_for(log_ixs[:, :, imode]),
                   facecolors='none', edgecolors='red', linewidths=1.4)

        mean_f = np.nanmean(freq[:, :, imode])
        ax.set_title(f'Mode {imode + 1} (~{mean_f:.1f} {freq_label})')
        ax.set_xlabel(r'$\Delta$tth (deg)')
        ax.set_ylabel(r'$\Delta\gamma$ (deg)')

    for k in range(n_modes, nrows * ncols):
        fig.delaxes(axes[k // ncols, k % ncols])

    # Marker-size legend for log10(IXS_stokes)
    legend_vals = np.linspace(vmin, vmax, 3)
    legend_handles = [
        Line2D([], [], marker='o', linestyle='', markerfacecolor='none',
               markeredgecolor='red', markeredgewidth=1.4,
               markersize=np.sqrt(sizes_for(np.array([v]))[0]),
               label=f'$10^{{{v:.1f}}}$')
        for v in legend_vals
    ]
    fig.legend(handles=legend_handles, ncol=len(legend_handles),
               title='IXS Stokes (barn/sr, circle area ~ log scale)',
               loc='lower center', frameon=False)

    title = ('Analyzer Array — Phonon Frequencies (contours) '
             '& IXS Stokes Intensity (circles)')
    if Q_label:
        title += f'\nQ_center = {Q_label}'
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    return fig
