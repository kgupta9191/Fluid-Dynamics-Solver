import sys
import types
from pathlib import Path

_FAKE_MPI_INSTALLED = False


try:
    from mpi4py import MPI  # noqa: F401
except (ImportError, RuntimeError):
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
    _FAKE_MPI_INSTALLED = True


def pytest_sessionfinish(session, exitstatus):
    generated_file = Path.cwd() / "pressureSolverConvergenceData.txt"
    if generated_file.exists():
        generated_file.unlink()
    if _FAKE_MPI_INSTALLED:
        sys.modules.pop("mpi4py", None)
