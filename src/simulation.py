from input import *
from bc import *
from interVel import *
from press import *
from finalVel import *
from mpi4py import MPI
import constant
import time
import h5py

rank = comm.Get_rank()
size = comm.Get_size()
t = constant.Start_time # Start time for simulation
realTime = constant.End_time
iter = realTime/del_t
tic = time.time()
plt = 0
pldt = 1
while t < iter:
    u, v = boundaryCondition(u, v, u_bc, v_bc, u_t, u_b, u_l, u_r, v_t, v_b, v_l, v_r, nx, ny) # Set boundary conditions
    u, v = intermediateVelocity(u, v, rho, mu, nx, ny, dx, dy, del_t) #Solve for intermediate velocity condition
    p = pressureSolve(u, v, rho, nx, ny, dx, dy, del_t, iter) #Pressure iterative solver
    u, v, u_bc, v_bc = finalVelocity(p, u, v, rho, nx, ny, dx, dy, del_t) #Advance time step
    
    
    
    plt = plt + 1
    t = t + 1
    realTime = realTime + del_t


toc = time.time()
dt = toc - tic

if rank == 0:
    print("Execution time for ",size, "cores:", dt, "s")
    # Open an HDF5 file for writing
    with h5py.File("output.cgns", "w") as f:
        # Create the required datasets
        f.create_dataset("Base", data=np.array([3]))
        f.create_dataset("Zone", data=np.array([nx, ny, 1]))
        f.create_dataset("Pressure", data=p)
        f.create_dataset("VelocityU", data=u)
        f.create_dataset("VelocityV", data=v)



U = np.zeros((ny, nx))
V = np.zeros((ny, nx))
velMag = np.zeros((ny, nx))

# Interpolate my cell centered U and V velocites
for i in range(nxc+1, nx):
    for j in range(1, ny):
        U[j, i] = 0.5 * (u[j+1, i+1] + u[j+1, i+2])
        V[j, i] = 0.5 * (v[j+1, i+1] + v[j+2, i+1])
        velMag[j, i] = np.sqrt(U[j, i]**2 + V[j, i]**2)
for i in range(1, nxc):
    for j in range(nyc+1, ny):
        U[j, i] = 0.5 * (u[j+1, i+1] + u[j+1, i+2])
        V[j, i] = 0.5 * (v[j+1, i+1] + v[j+2, i+1])
        velMag[j, i] = np.sqrt(U[j, i]**2 + V[j, i]**2)


if rank == 0:
    # Plot stream particles
    fig1, ax1 = plot.subplots(figsize=(6, 6))
    strm = ax1.streamplot(x, y, U, V, color=U, linewidth=2, cmap=plot.cm.autumn)
    cbar = plot.colorbar(strm.lines)
    cbar.ax.set_ylabel('Velocity magnitude')
    ax1.set_title('Stream Particles')
    ax1.set_xlabel('X-coordinate (m)')
    ax1.set_ylabel('Y-coordinate (m)')


    # Plot velocity vectors
    fig2, ax2 = plot.subplots()
    q = ax2.quiver(x, y, U, V, units='width', color='k', scale=2)
    ax2.set_title('Velocity Vector')
    ax2.set_xlabel('X-coordinate (m)')
    ax2.set_ylabel('Y-coordinate (m)')

plot.show()
