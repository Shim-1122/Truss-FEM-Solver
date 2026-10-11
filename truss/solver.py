"""Linear-elastic 2D truss solver using the direct stiffness method.

Sign convention: axial force > 0 is tension, < 0 is compression.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu

from .model import Truss


@dataclass
class Result:
    displacements: np.ndarray   # (n_nodes, 2) [m]
    axial_forces: np.ndarray    # (n_elements,) [N]
    stresses: np.ndarray        # (n_elements,) [Pa]
    reactions: np.ndarray       # (n_nodes, 2) support forces [N], zero at free DOFs


def element_geometry(truss: Truss) -> tuple[np.ndarray, np.ndarray]:
    """Length of every element, and its axial projection vector g = [-c, -s, c, s].

    (c, s) = (cos, sin) of the element angle. For element DOFs u_e = [ui, vi, uj, vj],
    the elongation is simply  g . u_e.
    """
    d = truss.nodes[truss.elements[:, 1]] - truss.nodes[truss.elements[:, 0]]
    length = np.hypot(d[:, 0], d[:, 1])
    c, s = d[:, 0] / length, d[:, 1] / length
    return length, np.column_stack([-c, -s, c, s])


def element_dofs(truss: Truss) -> np.ndarray:
    """(n_elements, 4) global DOF numbers of each element: [2i, 2i+1, 2j, 2j+1]."""
    i, j = truss.elements[:, 0], truss.elements[:, 1]
    return np.column_stack([2 * i, 2 * i + 1, 2 * j, 2 * j + 1])


def assemble_stiffness(truss: Truss) -> csr_matrix:
    """Global stiffness matrix K (sparse).

    Each element matrix is (EA/L) * g g^T (4x4). All are built in one vectorized
    step; the COO constructor sums duplicate (row, col) entries, which performs
    the assembly.
    """
    length, g = element_geometry(truss)
    k_axial = truss.E * truss.A / length                                  # (n_el,)
    k_el = k_axial[:, None, None] * g[:, :, None] * g[:, None, :]          # (n_el, 4, 4)

    dofs = element_dofs(truss)
    rows = np.broadcast_to(dofs[:, :, None], k_el.shape)
    cols = np.broadcast_to(dofs[:, None, :], k_el.shape)
    n = truss.n_dof
    return coo_matrix((k_el.ravel(), (rows.ravel(), cols.ravel())), shape=(n, n)).tocsr()


def solve(truss: Truss) -> Result:
    """Solve K u = F with zero displacement at the supports."""
    K = assemble_stiffness(truss)
    F = truss.loads.ravel()
    free = ~truss.fixed.ravel()

    # Partition K into free/fixed DOFs; supported DOFs are known (zero), so only K_ff u_f = F_f remains
    K_ff = K[free][:, free].tocsc()
    try:
        u_free = splu(K_ff).solve(F[free])
    except RuntimeError:                                    # exactly zero pivot
        u_free = None
    tolerance = 1e-6 * max(1.0, np.abs(F).max())
    if u_free is None or not np.allclose(K_ff @ u_free, F[free], atol=tolerance):
        raise ValueError("Stiffness matrix is singular: the truss is a mechanism. "
                         "Add supports or bracing members.")

    u = np.zeros(truss.n_dof)
    u[free] = u_free

    reactions = K @ u - F
    reactions[free] = 0.0                                   # only supports carry reactions

    length, g = element_geometry(truss)
    elongation = np.sum(g * u[element_dofs(truss)], axis=1)
    axial = truss.E * truss.A / length * elongation
    return Result(u.reshape(-1, 2), axial, axial / truss.A, reactions.reshape(-1, 2))
