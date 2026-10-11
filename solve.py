"""Solve a truss from a JSON file and report forces, reactions, and deflection.

    python solve.py examples/pratt_bridge.json
    python solve.py examples/pratt_bridge.json --scale 200 --out figures/bridge.png
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from truss import Truss, solve
from truss.plot import plot_truss


def report(truss: Truss, result) -> None:
    dof = truss.determinacy()
    status = "determinate" if dof == 0 else f"indeterminate, degree {dof}" if dof > 0 else "mechanism"
    print(f"{len(truss.nodes)} nodes, {len(truss.elements)} members, statically {status}\n")

    print(f"{'member':>6} {'nodes':>8} {'force (kN)':>11} {'stress (MPa)':>13}  state")
    for k, ((i, j), f, s) in enumerate(zip(truss.elements, result.axial_forces, result.stresses)):
        state = "tension" if f > 1e-6 else "compression" if f < -1e-6 else "zero-force"
        print(f"{k:>6} {f'{i}-{j}':>8} {f / 1e3:>11.2f} {s / 1e6:>13.2f}  {state}")

    print("\nReactions (kN):")
    for node in np.flatnonzero(truss.fixed.any(axis=1)):
        rx, ry = np.round(result.reactions[node] / 1e3, 2) + 0.0   # + 0.0 turns -0.00 into 0.00
        print(f"  node {node}: Rx = {rx:.2f}, Ry = {ry:.2f}")

    disp_mm = np.linalg.norm(result.displacements, axis=1) * 1e3
    worst = disp_mm.argmax()
    print(f"\nMax displacement: {disp_mm[worst]:.3f} mm at node {worst}")
    print(f"Max |stress|: {np.abs(result.stresses).max() / 1e6:.1f} MPa")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("file", help="truss definition (.json)")
    parser.add_argument("--scale", type=float, default=None, help="deformation magnification (default: auto)")
    parser.add_argument("--out", default=None, help="figure path (default: figures/<name>.png)")
    args = parser.parse_args()

    truss = Truss.from_json(args.file)
    result = solve(truss)
    report(truss, result)

    out = Path(args.out or f"figures/{Path(args.file).stem}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    name = json.loads(Path(args.file).read_text(encoding="utf-8")).get("name", Path(args.file).stem)
    ax = plot_truss(truss, result, args.scale, title=name)
    ax.figure.tight_layout()
    ax.figure.savefig(out, dpi=150)
    print(f"\nFigure saved to {out}")
