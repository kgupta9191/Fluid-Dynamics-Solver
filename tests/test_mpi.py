import numpy as np
from mpi4py import MPI

from input import del_t, dx, dy, nx, ny, rho
from press import pressureSolve


def test_mpi_world_properties():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    assert size >= 1
    assert 0 <= rank < size


def test_mpi_allreduce_sum_consistency():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    total = comm.allreduce(rank, op=MPI.SUM)
    expected = (size * (size - 1)) // 2
    assert total == expected


def test_pressure_solver_returns_finite_field_for_zero_velocity():
    u = np.zeros((ny + 2, nx + 2))
    v = np.zeros((ny + 2, nx + 2))

    p = pressureSolve(u, v, rho, nx, ny, dx, dy, del_t, t=0)

    assert p.shape == (ny + 2, nx + 2)
    assert np.isfinite(p).all()
