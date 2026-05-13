from input import *

def boundaryCondition(u, v, u_bc, v_bc, u_t, u_b, u_l, u_r, v_t, v_b, v_l, v_r, nx, ny):
    # Setting initial velocities (BC to be applied on ghost cells)
    # Bottom wall
    u_bc[0, :] = 2*u_b - u_bc[1, :]
    v_bc[1, :] = v_b
    # Left wall
    u_bc[:, 1] = u_l
    v_bc[:, 0] = 2*v_l - v_bc[:, 1] 
    # Right wall
    u_bc[:, nx+1] = u_r
    v_bc[:, nx+1] = 2*v_r - v_bc[:, nx]
    # Top wall
    u_bc[ny+1, :] = 2*u_t - u_bc[ny, :]
    v_bc[ny+1, :] = v_t
    # Cavity wall
    u_bc[0:nyc, 0:nxc+1] = u_c
    u_bc[nyc, 0:nxc] = 2*u_c - u_bc[nyc+1, 0:nxc]
    v_bc[0:nyc+1, 0:nxc] = v_c
    v_bc[0:nyc, nxc] = 2.0*v_c - v_bc[0:nyc, nxc+1]
    u[:] = u_bc[:]
    v[:] = v_bc[:]
    return u, v
