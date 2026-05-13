from input import *
import constant 
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def pressureSolve(u, v, rho, nx, ny, dx, dy, del_t, t):
    # Pressure loop
    # Initializing variables and pressure coefficients
    p = np.zeros((ny+2, nx+2))
    rhs = np.zeros((ny+2, nx+2))
    a_p = np.zeros((ny+2, nx+2))
    a_e = np.ones((ny+2, nx+2))/(rho*dx**2)
    a_w = np.ones((ny+2, nx+2))/(rho*dx**2)
    a_n = np.ones((ny+2, nx+2))/(rho*dy**2)
    a_s = np.ones((ny+2, nx+2))/(rho*dy**2)
    # Pressure coefficients along domain boundary
    a_e[:, nx] = 0.0
    a_w[:, 1] = 0.0
    a_n[ny, :] = 0.0
    a_s[1, :] = 0.0
    a_s[nyc+1, 1:nxc] = 0.0
    a_w[1:nyc, nxc+1] = 0.0
    # Central pressure coefficient
    a_p = -(a_e + a_w + a_n + a_s)
    # Pressure solver tolerance
    TOL = constant.Tolerance
    beta = 1.872
    counter = 0
    maxError = 10
    while (abs(maxError) > TOL):
        maxError = 0
        counter = counter + 1
        # Solve for pressure field
        chunk_size = nx / size
        start_i = int(rank * chunk_size) 
        end_i = int((rank + 1) * chunk_size)
        start_j, end_j = nyc, nx
        for i in range(start_i+1, end_i+1):
            for j in range(start_j+1, end_j+1):
                rhs[j, i] = +(u[j, i+1] - u[j, i])/dx + (v[j+1, i] - v[j, i])/dy
                rhs[j, i] = rhs[j, i]/del_t - (a_w[j, i]*p[j, i-1] + a_e[j, i]*p[j, i+1] + a_n[j, i]*p[j+1, i] + a_s[j, i]*p[j-1, i])
                p[j, i] = beta*rhs[j, i]/a_p[j, i] + (1-beta)*p[j, i]
        for i in range(start_i+1, end_i+1):
            for j in range(start_j+1, end_j+1):
                rhs[j, i] = +(u[j, i+1] - u[j, i])/dx + (v[j+1, i] - v[j, i])/dy
                error = a_w[j, i]*p[j, i-1] + a_e[j, i]*p[j, i+1] + a_n[j, i]*p[j+1, i] + a_s[j, i]*p[j-1, i] + a_p[j, i]*p[j, i] - (rhs[j, i]/del_t)
                if(np.abs(error) > maxError):
                    maxError = np.abs(error)
        chunk_size1 = (nx-nxc) / size
        start_ii = int(rank * chunk_size1) + nxc
        end_ii = int((rank + 1) * chunk_size1) + nxc
        start_jj, end_jj = 0, nyc
        for i in range(start_ii+1, end_ii+1):
            for j in range(start_jj+1, end_jj+1):
                rhs[j, i] = +(u[j, i+1] - u[j, i])/dx + (v[j+1, i] - v[j, i])/dy
                rhs[j, i] = rhs[j, i]/del_t - (a_w[j, i]*p[j, i-1] + a_e[j, i]*p[j, i+1] + a_n[j, i]*p[j+1, i] + a_s[j, i]*p[j-1, i])
                p[j, i] = beta*rhs[j, i]/a_p[j, i] + (1-beta)*p[j, i]
        for i in range(start_ii+1, end_ii+1):
            for j in range(start_jj+1, end_jj+1):
                rhs[j, i] = +(u[j, i+1] - u[j, i])/dx + (v[j+1, i] - v[j, i])/dy
                error = a_w[j, i]*p[j, i-1] + a_e[j, i]*p[j, i+1] + a_n[j, i]*p[j+1, i] + a_s[j, i]*p[j-1, i] + a_p[j, i]*p[j, i] - (rhs[j, i]/del_t)
                if(np.abs(error) > maxError):
                    maxError = np.abs(error)
        maxError = comm.allreduce(maxError, op=MPI.MAX)
    myFile.write('{:>8s} {:>8s} {:>8s}\n'.format(str(t), str(maxError), str(counter)))
    p = comm.allreduce(p)
    return p
