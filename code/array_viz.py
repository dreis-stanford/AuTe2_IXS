"""2D visualization of phonon frequencies and IXS intensity across the
BL43LXU analyzer array.

For each phonon branch, plots a contour map of the mode frequency over the
analyzer array's (delta_tth, delta_gamma) layout -- the angular offsets of
each analyzer from the array center, computed by
SixCircleInterface.analyzer_array_offsets(). The IXS Stokes cross-section at
each analyzer is overlaid as an open circle whose area encodes IXS_stokes on
a linear scale (saturating at both ends), so both quantities can be read off
the same plot.
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

    # IXS -> marker area, shared scale across all modes so sizes are
    # comparable subplot-to-subplot. IXS_stokes spans many decades (some
    # modes are selection-rule-forbidden, ~0), so a linear scale across the
    # full min/max would make almost everything either invisible or
    # saturated. Instead use a linear scale between the 10th/90th
    # percentiles of the (finite) values, saturating outside that range --
    # values at/below the low percentile floor out at the minimum marker
    # size ("a single point"), values at/above the high percentile cap out
    # at the maximum size.
    SATURATION_PERCENTILES = (10.0, 90.0)
    finite = ixs[np.isfinite(ixs)]
    if finite.size:
        vmin, vmax = np.percentile(finite, SATURATION_PERCENTILES)
    else:
        vmin, vmax = 0.0, 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0

    MIN_AREA = 3.0    # marker area, points^2
    MAX_AREA = 215.0  # marker area, points^2
    MIN_RADIUS = np.sqrt(MIN_AREA)
    MAX_RADIUS = np.sqrt(MAX_AREA)

    def sizes_for(vals):
        # Linear in radius (not area) so the visual size differences
        # between small and large IXS values are more pronounced.
        frac = np.clip((vals - vmin) / (vmax - vmin), 0.0, 1.0)
        radius = MIN_RADIUS + frac * (MAX_RADIUS - MIN_RADIUS)
        return radius**2  # marker area, points^2 (matplotlib scatter `s`)

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

        ax.scatter(dtth, dgam, s=sizes_for(ixs[:, :, imode]),
                   facecolors='none', edgecolors='red', linewidths=1.4)

        mean_f = np.nanmean(freq[:, :, imode])
        ax.set_title(f'Mode {imode + 1} (~{mean_f:.1f} {freq_label})')
        ax.set_xticks([])
        ax.set_yticks([])

    for k in range(n_modes, nrows * ncols):
        fig.delaxes(axes[k // ncols, k % ncols])

    # Marker-size legend for IXS_stokes (linear, saturating outside [vmin, vmax])
    legend_vals = np.linspace(vmin, vmax, 3)
    legend_labels = [f'≤{vmin:.2g}'] + [f'{v:.2g}' for v in legend_vals[1:-1]] \
        + [f'≥{vmax:.2g}']
    legend_handles = [
        Line2D([], [], marker='o', linestyle='', markerfacecolor='none',
               markeredgecolor='red', markeredgewidth=1.4,
               markersize=np.sqrt(sizes_for(np.array([v]))[0]),
               label=lbl)
        for v, lbl in zip(legend_vals, legend_labels)
    ]
    fig.legend(handles=legend_handles, ncol=len(legend_handles),
               title='IXS Stokes (barn/sr, circle area ~ linear scale, saturating)',
               loc='lower center', frameon=False)

    # Plot range, common to all subplots -- shown once here since the
    # per-subplot axes are unlabeled to save space.
    dtth_lo, dtth_hi = np.nanmin(dtth), np.nanmax(dtth)
    dgam_lo, dgam_hi = np.nanmin(dgam), np.nanmax(dgam)
    if grid is not None:
        dtth_lo, dtth_hi = grid_dtth.min(), grid_dtth.max()
        dgam_lo, dgam_hi = grid_dgam.min(), grid_dgam.max()

    title = ('Analyzer Array — Phonon Frequencies (contours) '
             '& IXS Stokes Intensity (circles)')
    if Q_label:
        title += f'\nQ_center = {Q_label}'
    title += (f'\nPlot range: '
              f'Δtth ∈ [{dtth_lo:.1f}, {dtth_hi:.1f}]°, '
              f'Δγ ∈ [{dgam_lo:.1f}, {dgam_hi:.1f}]°')
    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0.06, 1, 0.91))
    return fig
