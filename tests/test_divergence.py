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


def _divergence(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    return (u[1 : ny + 1, 2 : nx + 2] - u[1 : ny + 1, 1 : nx + 1]) / dx + (
        v[2 : ny + 2, 1 : nx + 1] - v[1 : ny + 1, 1 : nx + 1]
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
    p = pressureSolve(u_star.copy(), v_star.copy(), rho, nx, ny, dx, dy, del_t, t=0)
    u_next, v_next, _, _ = finalVelocity(
        p, u_star.copy(), v_star.copy(), rho, nx, ny, dx, dy, del_t
    )

    before = np.linalg.norm(_divergence(u_star, v_star))
    after = np.linalg.norm(_divergence(u_next, v_next))

    assert np.isfinite(before)
    assert np.isfinite(after)
    assert after <= before + 1e-10
