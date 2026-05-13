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


def test_single_step_solver_pipeline_runs_and_preserves_shapes():
    u = np.zeros((ny + 2, nx + 2))
    v = np.zeros((ny + 2, nx + 2))
    u_bc = np.zeros_like(u)
    v_bc = np.zeros_like(v)

    u, v = boundaryCondition(
        u, v, u_bc, v_bc, u_t, u_b, u_l, u_r, v_t, v_b, v_l, v_r, nx, ny
    )
    u, v = intermediateVelocity(u, v, rho, mu, nx, ny, dx, dy, del_t)
    p = pressureSolve(u, v, rho, nx, ny, dx, dy, del_t, t=0)
    u, v, u_bc, v_bc = finalVelocity(p, u, v, rho, nx, ny, dx, dy, del_t)

    assert u.shape == (ny + 2, nx + 2)
    assert v.shape == (ny + 2, nx + 2)
    assert p.shape == (ny + 2, nx + 2)
    assert u_bc.shape == (ny + 2, nx + 2)
    assert v_bc.shape == (ny + 2, nx + 2)

    assert np.isfinite(u).all()
    assert np.isfinite(v).all()
    assert np.isfinite(p).all()
