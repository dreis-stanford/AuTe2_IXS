"""3D visualization of the six-circle scattering geometry.

Draws the incoming/outgoing beams, the momentum transfer Q, the rotated
crystal axes, and the fixed lab coordinate system for a given set of
diffractometer angles, so you can see how the sample is oriented.

Conventions follow verify_scattering.py (the project's validated convention):
  x: horizontal, perpendicular to the beam
  y: beam axis (source at +y; the beam travels along -y to the sample at origin)
  z: vertical (up)
  k_in  = k [0, -cos(mu), -sin(mu)]                 (along -y, tilted by mu)
  k_out = k [sin(tth)cos(gam), -cos(tth)cos(gam), -sin(gam)]
  Q     = k_out - k_in
Crystal frame at all angles zero: a||x, b||-y, c||z; rotated into the lab frame
by the sample circles via R = R_mu @ R_theta @ R_chi @ R_phi.
"""
import numpy as np

from .verify_scattering import rotate_Q_to_lab_frame


def _kin_hat(mu_deg):
    mu = np.radians(mu_deg)
    return np.array([0.0, -np.cos(mu), -np.sin(mu)])


def _kout_hat(tth_deg, gam_deg):
    tth, gam = np.radians(tth_deg), np.radians(gam_deg)
    return np.array([np.sin(tth) * np.cos(gam),
                     -np.cos(tth) * np.cos(gam),
                     -np.sin(gam)])


def _arrow(ax, origin, vec, color, label=None, lw=2.5, alpha=1.0,
           label_scale=1.08):
    """Draw a 3D arrow from origin along vec, with an optional label at its tip."""
    o = np.asarray(origin, float)
    v = np.asarray(vec, float)
    ax.quiver(o[0], o[1], o[2], v[0], v[1], v[2],
              color=color, linewidth=lw, alpha=alpha,
              arrow_length_ratio=0.12)
    if label:
        tip = o + v * label_scale
        ax.text(tip[0], tip[1], tip[2], label, color=color,
                fontsize=11, fontweight='bold')


def _in_plane_basis(n, prefer):
    """Two orthonormal vectors spanning the plane perpendicular to unit n.
    The first follows `prefer` (projected into the plane) when possible."""
    n = n / np.linalg.norm(n)
    u = prefer - np.dot(prefer, n) * n
    if np.linalg.norm(u) < 1e-6:
        seed = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(seed, n)) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
        u = seed - np.dot(seed, n) * n
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    return u, v


def _slab_faces(center, u, v, n, hu, hv, ht):
    """Six quad faces of a rectangular slab centered at `center`, with
    half-extents hu, hv (in plane) and ht (along the thin normal n)."""
    c = np.asarray(center, float)
    def corner(su, sv, sn):
        return c + su * hu * u + sv * hv * v + sn * ht * n
    quads = []
    for sn in (-1, 1):  # broad faces (the mounted surfaces)
        quads.append([corner(-1, -1, sn), corner(1, -1, sn),
                      corner(1, 1, sn), corner(-1, 1, sn)])
    for su in (-1, 1):  # side faces
        quads.append([corner(su, -1, -1), corner(su, 1, -1),
                      corner(su, 1, 1), corner(su, -1, 1)])
    for sv in (-1, 1):
        quads.append([corner(-1, sv, -1), corner(1, sv, -1),
                      corner(1, sv, 1), corner(-1, sv, 1)])
    return quads


def _rodrigues(v, axis, angle_deg):
    """Rotate vector v about `axis` by angle_deg (right-handed)."""
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    t = np.radians(angle_deg)
    v = np.asarray(v, float)
    return (v * np.cos(t) + np.cross(axis, v) * np.sin(t)
            + axis * np.dot(axis, v) * (1 - np.cos(t)))


def _draw_angle_arc(ax, axis, ref, angle_deg, radius, color, label):
    """Draw an arc gauging a rotation of `angle_deg` about `axis`, starting
    from direction `ref`, at the given radius. Labels with the value."""
    if abs(angle_deg) < 0.5:
        return
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    ref = np.asarray(ref, float)
    ref = ref - np.dot(ref, axis) * axis    # project into the plane perp to axis
    if np.linalg.norm(ref) < 1e-9:
        return
    ref = ref / np.linalg.norm(ref) * radius
    ts = np.linspace(0, angle_deg, 40)
    pts = np.array([_rodrigues(ref, axis, t) for t in ts])
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=color, lw=1.6)
    ax.scatter(*pts[-1], color=color, s=14)          # marks the current angle
    mid = pts[len(pts) // 2] * 1.18
    ax.text(mid[0], mid[1], mid[2], label, color=color,
            fontsize=9, fontweight='bold')


def plot_scattering_geometry(angles, hkl=None, surface_normal=None, ub=None,
                             show_diffractometer=True, title=None,
                             save=None, show=True, ax=None):
    """Render the scattering geometry for a set of diffractometer angles.

    Parameters
    ----------
    angles : dict with keys tth, th, chi, phi, mu, gam (degrees)
    hkl    : optional (H, K, L) shown in the title
    title  : optional title override
    save   : optional path to write a PNG
    show   : call plt.show() (blocks until the window is closed)
    ax     : optional existing 3D axes to draw into

    Returns the matplotlib Axes.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3d proj)

    tth = float(angles['tth']); th = float(angles['th'])
    chi = float(angles['chi']); phi = float(angles['phi'])
    mu = float(angles.get('mu', 0.0)); gam = float(angles.get('gam', 0.0))

    # Wavevectors in units of k (=2pi/lambda). |k_in| = |k_out| = 1.
    k_in = _kin_hat(mu)
    k_out = _kout_hat(tth, gam)
    Q = k_out - k_in

    if surface_normal is None:
        try:
            from . import config
            surface_normal = config.SURFACE_NORMAL
        except Exception:
            surface_normal = (0, 0, 1)
    sn = np.array(surface_normal, float)

    def _to_lab(v):
        v = np.asarray(v, float)
        v = rotate_Q_to_lab_frame(v / np.linalg.norm(v), chi, phi, th, mu)
        return v / np.linalg.norm(v)

    if ub is not None:
        # Faithful frame from the actual UB (same orientation the angles use).
        # Columns of UB are a*,b*,c* (phi frame); inv(UB).T columns are the real
        # (generally non-orthogonal) lattice vectors a,b,c. The surface normal is
        # a reciprocal vector: c* = UB @ (h,k,l).
        ub = np.asarray(ub, float).reshape(3, 3)
        real_axes = np.linalg.inv(ub).T
        crystal = {nm: _to_lab(real_axes[:, i]) for i, nm in enumerate('abc')}
        n_lab = _to_lab(ub @ sn)               # SURFACE_NORMAL is in reciprocal space
        a_phi = real_axes[:, 0] / np.linalg.norm(real_axes[:, 0])
    else:
        # Fallback (no UB): idealized mounting a||+x, b||-y, c||+z. Approximate;
        # the crystal/normal will not match a real experimental UB.
        nominal = {'a': np.array([1.0, 0.0, 0.0]),
                   'b': np.array([0.0, -1.0, 0.0]),
                   'c': np.array([0.0, 0.0, 1.0])}
        crystal = {name: _to_lab(v) for name, v in nominal.items()}
        n_nom = (sn[0] * nominal['a'] + sn[1] * nominal['b'] + sn[2] * nominal['c'])
        n_lab = _to_lab(n_nom if np.linalg.norm(n_nom) > 1e-9 else nominal['c'])
        a_phi = nominal['a']
    # The surface normal is a FIXED mounting property (outward = +c*, the exposed
    # face) and is never flipped per reflection; flip the SURFACE_NORMAL HKL sign
    # in config if the opposite face is mounted out.

    if ax is None:
        fig = plt.figure(figsize=(9.0, 8.0))
        ax = fig.add_subplot(111, projection='3d')

    axis_colors = {'a': '#1b9e77', 'b': '#7570b3', 'c': '#d95f02'}
    c_kin, c_kout, c_Q, c_norm = '#e6194b', '#4363d8', 'black', '#e7298a'

    # --- Rectangular crystal slab, centered on the sample -------------------
    # Thin along the surface normal; broad faces follow the in-plane crystal
    # axes (a preferred for the long edge), so you see the crystal rotate.
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    u, v = _in_plane_basis(n_lab, prefer=crystal['a'])
    ht = 0.08
    # The exposed (mounted) face is the +n (+c*) face; offset the slab so it
    # passes through the sample (origin) and shade it. For (0,0,L>0) Q || +c*, so
    # this is also the reflecting face the beam strikes.
    slab_center = -ht * n_lab
    faces = _slab_faces(slab_center, u, v, n_lab, hu=0.5, hv=0.34, ht=ht)
    front_face = faces[1]                     # faces[1] = +n (+c*) broad face
    bulk_faces = [f for i, f in enumerate(faces) if i != 1]
    ax.add_collection3d(Poly3DCollection(
        bulk_faces, facecolor='#9ecae1', edgecolor='0.25',
        linewidths=0.8, alpha=0.30))
    ax.add_collection3d(Poly3DCollection(
        [front_face], facecolor='#fd8d3c', edgecolor='#b35806',
        linewidths=1.8, alpha=0.85))

    # Crystal axes (rotated) anchored at the sample.
    for name, vv in crystal.items():
        _arrow(ax, (0, 0, 0), 0.8 * vv, color=axis_colors[name],
               label=name, lw=3.0, label_scale=1.12)
    # Surface normal.
    _arrow(ax, (0, 0, 0), 0.95 * n_lab, color=c_norm, lw=2.2, label='n')

    # --- Beams and Q (identified in the legend to avoid label clutter) ------
    _arrow(ax, -k_in, k_in, color=c_kin, lw=2.5)        # incident, into sample
    _arrow(ax, (0, 0, 0), k_out, color=c_kout, lw=2.5)  # scattered, out of sample
    _arrow(ax, (0, 0, 0), Q, color=c_Q, lw=3.5)         # momentum transfer
    ax.plot([k_in[0], k_out[0]], [k_in[1], k_out[1]], [k_in[2], k_out[2]],
            color=c_Q, ls=':', lw=1.2, alpha=0.6)        # triangle closure
    ax.scatter([0], [0], [0], color='black', s=30)

    # --- Simplified diffractometer (schematic, not to scale) ----------------
    # Shows how chi, theta rotate; tth is read from k_out (no tth circle), and a
    # small detector pad sits on the k_out arm. Axes/arcs are thin and live in
    # different planes to keep the figure uncluttered.
    c_dial = '#555555'
    if show_diffractometer:
        # Shared vertical axis for theta (sample) and 2theta (detector arm).
        ax.plot([0, 0], [0, 0], [-1.1, 1.1], color=c_dial, ls='--', lw=1.0)
        # theta: sample rotation about +z, gauged from the beam direction (-y).
        _draw_angle_arc(ax, axis=[0, 0, 1], ref=[0, -1, 0], angle_deg=th,
                        radius=0.55, color='#1f78b4', label=f'θ={th:.1f}°')
        # chi: sample tilt about the beam axis (-y), gauged from vertical (+z).
        _draw_angle_arc(ax, axis=[0, -1, 0], ref=[0, 0, 1], angle_deg=chi,
                        radius=0.8, color='#33a02c', label=f'χ={chi:.1f}°')
        # phi: innermost sample rotation, about the chi/theta-tilted sample axis
        # (z carried by mu,theta,chi), gauged from the phi=0 position of a.
        phi_axis = rotate_Q_to_lab_frame([0, 0, 1], chi, 0.0, th, mu)
        phi_ref = rotate_Q_to_lab_frame(a_phi, chi, 0.0, th, mu)
        _draw_angle_arc(ax, axis=phi_axis, ref=phi_ref, angle_deg=phi,
                        radius=0.98, color='#b15928', label=f'φ={phi:.1f}°')
        # 2theta is read directly from the k_out direction (no tth circle/pad).

    # --- Fixed lab frame, offset to the (empty) upper-right corner ----------
    lim = 1.4
    lab_origin = np.array([lim * 0.85, lim * 0.85, lim * 0.80])
    lab_len = 0.35  # deliberately smaller than k_in (=1): just a reference
    for axis_vec, name in [((lab_len, 0, 0), 'x'), ((0, lab_len, 0), 'y'),
                           ((0, 0, lab_len), 'z')]:
        _arrow(ax, lab_origin, axis_vec, color='0.25', lw=1.8,
               label=name, label_scale=1.25)

    # Legend (high-contrast on white).
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], color=c_kin, lw=2.5, label='k_in (incident)'),
        Line2D([0], [0], color=c_kout, lw=2.5, label='k_out (scattered)'),
        Line2D([0], [0], color=c_Q, lw=3.5, label='Q = k_out - k_in'),
        Line2D([0], [0], color=axis_colors['a'], lw=3.0, label='crystal a'),
        Line2D([0], [0], color=axis_colors['b'], lw=3.0, label='crystal b'),
        Line2D([0], [0], color=axis_colors['c'], lw=3.0, label='crystal c'),
        Line2D([0], [0], color=c_norm, lw=2.2,
               label=f'surface normal ({surface_normal[0]:g}{surface_normal[1]:g}{surface_normal[2]:g})'),
        Patch(facecolor='#9ecae1', edgecolor='0.25', alpha=0.5, label='crystal'),
        Patch(facecolor='#fd8d3c', edgecolor='0.25', alpha=0.5, label='front face (beam)'),
        Line2D([0], [0], color='0.25', lw=1.8, label='lab x, y, z (corner)'),
    ]
    ax.legend(handles=handles, loc='upper left', fontsize=9,
              framealpha=0.9, edgecolor='0.3')

    # Cosmetics: equal aspect, fixed limits, view onto the scattering plane.
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
    ax.view_init(elev=22, azim=-66)

    if title is None:
        hkl_str = f"  HKL=({hkl[0]:g}, {hkl[1]:g}, {hkl[2]:g})" if hkl is not None else ""
        title = ("Scattering geometry" + hkl_str + "\n"
                 f"tth={tth:.2f}  th={th:.2f}  chi={chi:.2f}  "
                 f"phi={phi:.2f}  mu={mu:.3f}  gam={gam:.2f}   "
                 f"|Q|={np.linalg.norm(Q):.3f} k")
    ax.set_title(title, fontsize=10)

    if save:
        import os
        d = os.path.dirname(save)
        if d:
            os.makedirs(d, exist_ok=True)
        ax.figure.savefig(save, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    return ax


def _main(argv=None):
    """CLI entry point so the plot can run in its own process (interactive,
    rotatable 3D window) without blocking a calling REPL.

    Example:
        python -m code.geometry_viz --angles 7.89,3.94,55.08,-0.02,-0.17,0 \\
                                    --hkl 1,0,1
    """
    import argparse
    p = argparse.ArgumentParser(description="3D scattering geometry view")
    p.add_argument('--angles', required=True,
                   help='tth,th,chi,phi,mu,gam in degrees (comma separated)')
    p.add_argument('--hkl', default=None, help='H,K,L (comma separated)')
    p.add_argument('--surface-normal', default=None, help='H,K,L of surface normal')
    p.add_argument('--ub', default=None, help='UB matrix, 9 comma-separated floats (row-major)')
    p.add_argument('--save', default=None, help='optional PNG path')
    args = p.parse_args(argv)

    keys = ['tth', 'th', 'chi', 'phi', 'mu', 'gam']
    angles = dict(zip(keys, (float(x) for x in args.angles.split(','))))
    hkl = tuple(float(x) for x in args.hkl.split(',')) if args.hkl else None
    sn = tuple(float(x) for x in args.surface_normal.split(',')) \
        if args.surface_normal else None
    ub = [float(x) for x in args.ub.split(',')] if args.ub else None
    plot_scattering_geometry(angles, hkl=hkl, surface_normal=sn, ub=ub,
                             save=args.save, show=True)


if __name__ == '__main__':
    _main()
