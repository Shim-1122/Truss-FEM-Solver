"""Draw a truss: undeformed shape, exaggerated deformed shape, members colored by axial force."""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from .model import Truss
from .solver import Result


def plot_truss(truss: Truss, result: Result, scale: float | None = None, title: str = "", ax=None, label_forces=True):
    """Red = tension, blue = compression. `scale` magnifies displacements (auto if None)."""
    if ax is None:
        _, ax = plt.subplots(figsize=(11, 5.5))

    size = np.ptp(truss.nodes, axis=0).max()
    max_disp = np.abs(result.displacements).max()
    if scale is None:
        raw = 0.03 * size / max_disp if max_disp > 0 else 1.0     # largest displacement drawn at 3% of size
        scale = float(f"{raw:.1g}")                                 # round to 1 significant figure
    deformed = truss.nodes + scale * result.displacements

    # Undeformed shape in light gray, deformed shape colored by force
    segments = lambda xy: xy[truss.elements]                        # (n_el, 2 endpoints, 2 coords)
    ax.add_collection(LineCollection(segments(truss.nodes), colors="0.8", linewidths=1.5, linestyles="--"))
    kN = result.axial_forces / 1e3
    limit = np.abs(kN).max() or 1.0
    members = LineCollection(segments(deformed), array=kN, cmap="coolwarm",
                             norm=Normalize(-limit, limit), linewidths=3.5)
    ax.add_collection(members)
    plt.colorbar(members, ax=ax, orientation="horizontal", fraction=0.06, pad=0.18, aspect=40,
                 label="axial force (kN), tension positive")

    if label_forces:
        for (p, q), f in zip(deformed[truss.elements], kN):
            ax.text(*(p + q) / 2, f"{f:.0f}" if abs(f) >= 0.5 else "0",   # avoid printing "-0" fontsize=7, ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))

    ax.plot(*deformed.T, "o", color="k", ms=4, zorder=3)
    for node in np.flatnonzero(truss.fixed.any(axis=1)):             # supports
        marker = "^" if truss.fixed[node].all() else "o"
        ax.plot(*truss.nodes[node], marker, ms=14, mfc="none", mec="k", mew=1.5, zorder=2)
    arrow_len = 0.08 * size
    for node in np.flatnonzero(np.abs(truss.loads).sum(axis=1) > 0):  # applied loads
        direction = truss.loads[node] / np.linalg.norm(truss.loads[node])
        tail = deformed[node] - direction * arrow_len
        ax.annotate("", xy=deformed[node], xytext=tail, arrowprops=dict(arrowstyle="-|>", color="C2", lw=2))

    ax.set_aspect("equal")
    ax.autoscale()
    ax.margins(0.08)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"{title} (Deformation x{scale:g})")
    return ax
