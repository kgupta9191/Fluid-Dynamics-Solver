# Fluid Dynamics Solver (Incompressible Navier–Stokes)

A Python-based 2D incompressible CFD solver for a lid-driven cavity style benchmark, with pressure-projection time marching, MPI-enabled pressure solve, and visualization/output support.

---

## Table of Contents

- [Overview](#overview)
- [Physics and Numerical Method](#physics-and-numerical-method)
- [Code Architecture](#code-architecture)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Quick Start](#quick-start)
- [Running with MPI](#running-with-mpi)
- [Configuration](#configuration)
- [Outputs](#outputs)
- [Testing](#testing)
- [Current Limitations](#current-limitations)
- [Roadmap Ideas](#roadmap-ideas)
- [License](#license)

---

## Overview

This project solves the 2D incompressible Navier–Stokes equations on a structured Cartesian grid for a cavity-flow setup driven by a moving lid.

At a high level, each timestep performs:
1. Boundary condition update
2. Intermediate velocity update (momentum equations)
3. Pressure Poisson solve (projection step)
4. Final divergence-reduced velocity correction

The solver is written in plain NumPy and organized as modular stages, making it approachable for learning, research prototypes, and algorithm experimentation.

---

## Physics and Numerical Method

### Governing equations

The solver targets incompressible flow:

- Continuity: `∇·u = 0`
- Momentum: `∂u/∂t + (u·∇)u = -(1/ρ)∇p + ν∇²u`

where `u=(u,v)`, `ρ` is density, `ν=μ/ρ`, and `p` is pressure.

### Discretization and algorithm

- Structured 2D grid with ghost cells
- Explicit update for convection/diffusion in momentum predictor step
- Pressure correction via iterative Poisson solve
- Velocity projection to reduce divergence
- Boundary conditions enforced every timestep

### Boundary conditions

Default setup is a lid-driven cavity style flow:

- Top wall: moving tangential velocity (`u_t = Initial_Velocity`, `v_t = 0`)
- Other walls: no-slip / no-penetration style conditions
- Additional cavity/cutout boundary treatment is included in `bc.py` using `Cavity_cells_in_x_direction` and `Cavity_cells_in_y_direction`

---

## Code Architecture

Main pipeline (`src/simulation.py`):

- `boundaryCondition(...)` in `src/bc.py`
  - Applies wall and cavity ghost-cell boundary updates
- `intermediateVelocity(...)` in `src/interVel.py`
  - Computes momentum predictor (`u*`, `v*`)
- `pressureSolve(...)` in `src/press.py`
  - Iterative pressure solve with MPI all-reduce synchronization
- `finalVelocity(...)` in `src/finalVel.py`
  - Applies pressure gradient correction to get next-step velocities

Global setup and shared fields live in `src/input.py` and `src/constant.py`.

---

## Repository Structure

```text
Fluid-Dynamics-Solver/
├── src/
│   ├── simulation.py      # Main time loop and plotting
│   ├── input.py           # Grid, fields, timestep, and BC setup
│   ├── constant.py        # User-editable physical/numerical constants
│   ├── bc.py              # Boundary condition logic
│   ├── interVel.py        # Intermediate velocity (predictor) step
│   ├── press.py           # Pressure solver (MPI-aware)
│   └── finalVel.py        # Velocity correction step
├── tests/                 # Pytest-based validation suite
├── reports/report.pdf     # Project report
├── requirements.txt       # Python dependencies
├── run.sh                 # Convenience runner for local execution
└── README.md
```

---

## Requirements

From `requirements.txt`:

- numpy
- matplotlib
- mpi4py
- h5py
- pytest

Also required:

- Python 3.10+ (CI uses 3.10)
- MPI runtime for parallel runs (`mpiexec` / OpenMPI or equivalent)

---

## Quick Start

### 1) Install dependencies

```bash
pip install -r requirements.txt
```

### 2) Run the solver directly

```bash
python src/simulation.py
```

### 3) Or use the helper script

```bash
bash run.sh
```

`run.sh` checks for Python and dependencies, then runs `src/simulation.py`.

---

## Running with MPI

The pressure solve uses MPI collectives, so you can launch with multiple ranks:

```bash
mpiexec -n 4 python src/simulation.py
```

After completion, rank 0 prints timing and writes output.

---

## Configuration

Edit `src/constant.py` to control simulation setup.

### Key parameters

| Parameter | Meaning | Default |
|---|---|---|
| `Geometry_in_x_direction` | Domain length in x | `1.0` |
| `Geometry_in_y_direction` | Domain length in y | `1.0` |
| `Cells_in_x_direction` | Grid cells in x | `12` |
| `Cells_in_y_direction` | Grid cells in y | `12` |
| `Cavity_cells_in_x_direction` | Cavity/cutout size in x | `2` |
| `Cavity_cells_in_y_direction` | Cavity/cutout size in y | `2` |
| `Density` | Fluid density `ρ` | `1.0` |
| `Viscosity` | Dynamic viscosity `μ` | `0.01` |
| `Initial_Velocity` | Lid velocity | `1.0` |
| `Tolerance` | Pressure-solver convergence tolerance | `1e-2` |
| `Start_time` | Start timestep index | `0` |
| `End_time` | End timestep index | `1` |

`src/input.py` computes derived values (`dx`, `dy`, `del_t`) from these constants.

---

## Outputs

On run completion, the solver produces:

- `output.cgns` (HDF5-backed file written via `h5py`) containing:
  - `Pressure`
  - `VelocityU`
  - `VelocityV`
  - basic metadata datasets (`Base`, `Zone`)
- `pressureSolverConvergenceData.txt` with per-step residual and iteration data
- Matplotlib visualizations:
  - streamplot
  - quiver (velocity vectors)

---

## Testing

Run the test suite from repository root:

```bash
python -m pytest -v
```

Current tests cover:

- solver pipeline smoke behavior and field-shape checks
- projection-step divergence reduction behavior
- MPI world/allreduce sanity
- pressure solve finite-value behavior for zero-velocity input

---

## Current Limitations

- Baseline implementation is educational/research-oriented (not performance-optimized)
- Pressure solver and loop structure can be further optimized/vectorized
- Plotting is interactive (`matplotlib`) and may need adaptation for headless environments
- Output file uses HDF5 datasets in a CGNS-like container name, not a full validated CGNS hierarchy/toolchain workflow

---

## Roadmap Ideas

- Config file support (YAML/JSON/TOML) to avoid source edits for parameter sweeps
- Structured benchmark cases (e.g., Re=100, 400, 1000) and reference comparisons
- Higher-order advection options and stability/performance studies
- Improved domain decomposition and communication pattern for pressure solve
- Post-processing notebooks and automated report generation

---

## License

This repository is licensed under the terms in [LICENSE](LICENSE).
