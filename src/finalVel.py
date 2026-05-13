from input import *
from press import *
from bc import *

def finalVelocity(p, u, v, rho, nx, ny, dx, dy, del_t):
    for i in range (nxc+2, nx+1):
        for j in range (1, ny+1):
            # Update U velocity
            u[j, i] -= del_t*(p[j, i] - p[j, i-1])/(dx*rho)
    for i in range (2, nxc+2):
        for j in range (nyc+1, ny+1):
            # Update U velocity
            u[j, i] -= del_t*(p[j, i] - p[j, i-1])/(dx*rho)
    for i in range (nxc+1, nx+1):
        for j in range (2, ny+1):
            # Update V velocity
            v[j, i] -= del_t*(p[j, i] - p[j-1, i])/(dy*rho)
    for i in range (1, nxc+1):
        for j in range (nyc+2, ny+1):
            # Update V velocity
            v[j, i] -= del_t*(p[j, i] - p[j-1, i])/(dy*rho)
    u_bc[:] = u[:]
    v_bc[:] = v[:]
    return u, v, u_bc, v_bc
