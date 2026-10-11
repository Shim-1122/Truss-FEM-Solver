"""Truss definition: geometry, material, supports, and loads."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Support types -> which directions are fixed (x, y)
SUPPORT_TYPES = {"pin": (True, True), "roller": (False, True), "roller_x": (True, False)}


@dataclass
class Truss:
    nodes: np.ndarray       # (n_nodes, 2) coordinates [m]
    elements: np.ndarray    # (n_elements, 2) node index pairs
    E: np.ndarray           # (n_elements,) Young's modulus [Pa]
    A: np.ndarray           # (n_elements,) cross-section area [m^2]
    fixed: np.ndarray       # (n_nodes, 2) bool, True where displacement is prevented (x, y)
    loads: np.ndarray       # (n_nodes, 2) applied nodal force [N]

    @classmethod
    def build(cls, nodes, elements, E, A, supports: dict, loads: dict | None = None) -> "Truss":
        """Build from plain Python data.

        supports : {node: "pin" | "roller" | "roller_x" | [fix_x, fix_y]}
        loads    : {node: [Fx, Fy]}
        E, A     : one value for every member, or a list with one value per member
        """
        nodes = np.asarray(nodes, dtype=float)
        elements = np.asarray(elements, dtype=int)
        n_el = len(elements)

        fixed = np.zeros((len(nodes), 2), dtype=bool)
        for node, kind in supports.items():
            fixed[int(node)] = SUPPORT_TYPES[kind] if isinstance(kind, str) else kind

        force = np.zeros((len(nodes), 2))
        for node, f in (loads or {}).items():
            force[int(node)] = f

        E = np.broadcast_to(np.asarray(E, dtype=float), n_el).copy()
        A = np.broadcast_to(np.asarray(A, dtype=float), n_el).copy()
        return cls(nodes, elements, E, A, fixed, force)

    @classmethod
    def from_json(cls, path) -> "Truss":
        spec = json.loads(Path(path).read_text(encoding="utf-8"))
        return cls.build(spec["nodes"], spec["elements"], spec["E"], spec["A"],
                         spec["supports"], spec.get("loads"))

    @property
    def n_dof(self) -> int:
        return self.nodes.size          # 2 degrees of freedom (x, y) per node

    def determinacy(self) -> int:
        """Degree of static indeterminacy, m + r - 2j.
        0 = determinate, > 0 = indeterminate (needs stiffness, like this solver), < 0 = mechanism."""
        return len(self.elements) + int(self.fixed.sum()) - self.n_dof
