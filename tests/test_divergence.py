import numpy as np

from bc import boundaryCondition
from finalVel import finalVelocity
from input import (
    del_t,
    dx,
    dy,
    mu,
    nx,
    ny,
    rho,
    u_b,
    u_l,
    u_r,
    u_t,
    v_b,
    v_l,
    v_r,
    v_t,
)
from interVel import intermediateVelocity
from press import pressureSolve

DIVERGENCE_TOLERANCE = 1e-10

def _divergence(u: np.ndarray, v: np.ndarray, grid_nx: int, grid_ny: int) -> np.ndarray:
    return (u[1 : grid_ny + 1, 2 : grid_nx + 2] - u[1 : grid_ny + 1, 1 : grid_nx + 1]) / dx + (
        v[2 : grid_ny + 2, 1 : grid_nx + 1] - v[1 : grid_ny + 1, 1 : grid_nx + 1]
    ) / dy


def test_projection_step_reduces_divergence():
    u = np.zeros((ny + 2, nx + 2))
    v = np.zeros((ny + 2, nx + 2))
    u_bc = np.zeros_like(u)
    v_bc = np.zeros_like(v)

    u, v = boundaryCondition(
        u, v, u_bc, v_bc, u_t, u_b, u_l, u_r, v_t, v_b, v_l, v_r, nx, ny
    )
    u_star, v_star = intermediateVelocity(u.copy(), v.copy(), rho, mu, nx, ny, dx, dy, del_t)
    before = np.linalg.norm(_divergence(u_star, v_star, nx, ny))

    p = pressureSolve(u_star, v_star, rho, nx, ny, dx, dy, del_t, t=0)
    u_next, v_next, _, _ = finalVelocity(p, u_star, v_star, rho, nx, ny, dx, dy, del_t)
    after = np.linalg.norm(_divergence(u_next, v_next, nx, ny))

    assert np.isfinite(before)
    assert np.isfinite(after)
    assert after <= before + DIVERGENCE_TOLERANCE
