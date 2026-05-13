import sys
import types
from pathlib import Path


try:
    from mpi4py import MPI  # noqa: F401
except Exception:
    class _FakeComm:
        def Get_size(self):
            return 1

        def Get_rank(self):
            return 0

        def allreduce(self, value, op=None):
            return value

    class _FakeMPI:
        COMM_WORLD = _FakeComm()
        SUM = "sum"
        MAX = "max"

    fake_mpi4py = types.ModuleType("mpi4py")
    fake_mpi4py.MPI = _FakeMPI()
    sys.modules["mpi4py"] = fake_mpi4py


def pytest_sessionfinish(session, exitstatus):
    generated_file = Path.cwd() / "pressureSolverConvergenceData.txt"
    if generated_file.exists():
        generated_file.unlink()
