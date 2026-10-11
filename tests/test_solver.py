"""Solver checks against hand calculations and closed-form mechanics."""
from pathlib import Path

import numpy as np
import pytest

from truss import Truss, assemble_stiffness, solve

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
E, A = 200e9, 1e-4


def member(truss, i, j):
    """Index of the member connecting nodes i and j (either order)."""
    return next(k for k, (a, b) in enumerate(truss.elements) if {a, b} == {i, j})


# ------------------------------------------------------------ single members
def test_single_bar_elongation_is_PL_over_EA():
    P, L = 10e3, 2.0
    t = Truss.build([[0, 0], [L, 0]], [[0, 1]], E, A, {0: "pin", 1: "roller"}, {1: [P, 0]})
    r = solve(t)
    assert r.displacements[1, 0] == pytest.approx(P * L / (E * A))
    assert r.axial_forces[0] == pytest.approx(P)


@pytest.mark.parametrize("angle_deg", [0, 30, 90, 137, 250])
def test_result_does_not_depend_on_orientation(angle_deg):
    """Rotating the whole problem must not change the member force."""
    c, s = np.cos(np.radians(angle_deg)), np.sin(np.radians(angle_deg))
    tip = 3.0 * np.array([c, s])
    t = Truss.build([[0, 0], tip, tip + [-s, c]], [[0, 1], [1, 2], [0, 2]], E, A,
                    {0: "pin", 2: "pin"}, {1: list(5e3 * np.array([c, s]))})
    assert solve(t).axial_forces[0] == pytest.approx(5e3, rel=1e-9)


# --------------------------------------------------------- determinate truss
def test_two_bar_hanger_forces():
    """Load P hanging from two 45 degree bars: each carries P / (2 sin 45) in tension."""
    P = 1e3
    t = Truss.build([[0, 0], [2, 0], [1, -1]], [[0, 2], [1, 2]], E, A,
                    {0: "pin", 1: "pin"}, {2: [0, -P]})
    np.testing.assert_allclose(solve(t).axial_forces, P / (2 * np.sin(np.pi / 4)))


def test_pratt_bridge_matches_method_of_sections():
    t = Truss.from_json(EXAMPLES / "pratt_bridge.json")
    r = solve(t)
    # Reactions: 5 loads x 50 kN, symmetric -> 125 kN each
    assert r.reactions[0, 1] == pytest.approx(125e3)
    assert r.reactions[6, 1] == pytest.approx(125e3)
    # Section through panel 8-12 m, moments about top node 8 (x = 8): F = M(8)/h = (125*8 - 50*4)/4
    assert r.axial_forces[member(t, 2, 3)] == pytest.approx(200e3)
    # Moments about bottom node 3 (x = 12): F = -M(12)/h = -(125*12 - 50*4 - 50*8)/4
    assert r.axial_forces[member(t, 8, 9)] == pytest.approx(-225e3)
    # Unloaded center vertical is a zero-force member
    assert r.axial_forces[member(t, 3, 9)] == pytest.approx(0, abs=1e-6)


# ------------------------------------------------------- indeterminate truss
def test_three_bar_indeterminate_truss():
    """Classic statically indeterminate problem (needs compatibility, not just statics).
    Vertical bar: P / (1 + 2 cos^3 a).  Inclined bars: P cos^2 a / (1 + 2 cos^3 a)."""
    P, h, a = 10e3, 1.0, np.radians(40)
    x = h * np.tan(a)
    t = Truss.build([[-x, h], [0, h], [x, h], [0, 0]], [[0, 3], [1, 3], [2, 3]], E, A,
                    {0: "pin", 1: "pin", 2: "pin"}, {3: [0, -P]})
    assert t.determinacy() == 1
    f = solve(t).axial_forces
    denom = 1 + 2 * np.cos(a) ** 3
    assert f[1] == pytest.approx(P / denom)
    assert f[0] == pytest.approx(P * np.cos(a) ** 2 / denom)
    assert f[2] == pytest.approx(f[0])


def test_global_equilibrium_of_indeterminate_cantilever():
    t = Truss.from_json(EXAMPLES / "cantilever_x_braced.json")
    assert t.determinacy() == 4
    r = solve(t)
    total = r.reactions + t.loads
    np.testing.assert_allclose(total.sum(axis=0), 0, atol=1e-6)                    # sum F = 0
    moment = np.sum(t.nodes[:, 0] * total[:, 1] - t.nodes[:, 1] * total[:, 0])
    assert moment == pytest.approx(0, abs=1e-6)                                     # sum M = 0


# -------------------------------------------------------- stiffness matrix
def test_free_truss_has_three_rigid_body_modes():
    """Unsupported 2D body: 2 translations + 1 rotation cost zero energy -> 3 zero eigenvalues."""
    t = Truss.from_json(EXAMPLES / "pratt_bridge.json")
    K = assemble_stiffness(t).toarray()
    np.testing.assert_allclose(K, K.T)
    eig = np.linalg.eigvalsh(K)
    assert np.sum(np.abs(eig) < 1e-8 * eig.max()) == 3


def test_mechanism_is_detected():
    square = Truss.build([[0, 0], [1, 0], [1, 1], [0, 1]], [[0, 1], [1, 2], [2, 3], [3, 0]], E, A,
                         {0: "pin", 1: "roller"}, {2: [1e3, 0]})
    with pytest.raises(ValueError, match="mechanism"):
        solve(square)
