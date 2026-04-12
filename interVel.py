from input import *

def intermediateVelocity(u, v, rho, mu, nx, ny, dx, dy, del_t):
    # Solving N-S equation for momentum control volume in x-axis  
    for i in range (nxc+2, nx+1):
        for j in range (1, ny+1):
            u_e = 0.5*(u[j, i] + u[j, i+1])
            u_w = 0.5*(u[j, i-1] + u[j, i])
            u_n = 0.5*(u[j, i] + u[j+1, i])
            u_s = 0.5*(u[j, i] + u[j-1, i])
            v_n = 0.5*(v[j+1, i-1] + v[j+1, i])
            v_s = 0.5*(v[j, i-1] + v[j, i])
            convection = -((rho*u_e**2 - rho*u_w**2)/dx)-((rho*u_n*v_n - rho*u_s*u_s)/dy)
            diffusion = mu*((u[j, i+1]-2*u[j, i]+u[j, i-1])/dx**2 + (u[j+1, i]-2*u[j, i]+u[j-1, i])/dy**2)
            # calculating the intermediate velocity
            u[j, i] += (del_t*(convection+diffusion)/rho)  
    for i in range (2, nxc+2):
        for j in range (nyc+1, ny+1):
            u_e = 0.5*(u[j, i] + u[j, i+1])
            u_w = 0.5*(u[j, i-1] + u[j, i])
            u_n = 0.5*(u[j, i] + u[j+1, i])
            u_s = 0.5*(u[j, i] + u[j-1, i])
            v_n = 0.5*(v[j+1, i-1] + v[j+1, i])
            v_s = 0.5*(v[j, i-1] + v[j, i])
            convection = -((rho*u_e**2 - rho*u_w**2)/dx)-((rho*u_n*v_n - rho*u_s*u_s)/dy)
            diffusion = mu*((u[j, i+1]-2*u[j, i]+u[j, i-1])/dx**2 + (u[j+1, i]-2*u[j, i]+u[j-1, i])/dy**2)
            # calculating the intermediate velocity
            u[j, i] += (del_t*(convection+diffusion)/rho)  
    # Solving N-S equation for momentum control volume in y-axis 
    for i in range (nxc+1, nx+1):
        for j in range (2, ny+1):
            v_e = 0.5*(v[j, i] + v[j, i+1])
            v_w = 0.5*(v[j, i] + v[j, i-1])
            v_n = 0.5*(v[j, i] + v[j+1, i])
            v_s = 0.5*(v[j, i] + v[j-1, i])
            u_e = 0.5*(u[j-1, i+1] + u[j, i+1])
            u_w = 0.5*(u[j-1, i] + u[j, i])
            convection = -((rho*u_e*v_e - rho*u_w*v_w)/dx)-((rho*v_n**2 - rho*v_s**2)/dy)
            diffusion = mu*((v[j, i+1]-2*v[j, i]+v[j, i-1])/dx**2 + (v[j+1, i]-2*v[j, i]+v[j-1, i])/dy**2)
            # calculating the intermediate velocity
            v[j, i] += (del_t*(convection+diffusion)/rho) 
    for i in range (1, nxc+1):
        for j in range (nyc+2, ny+1):
            v_e = 0.5*(v[j, i] + v[j, i+1])
            v_w = 0.5*(v[j, i] + v[j, i-1])
            v_n = 0.5*(v[j, i] + v[j+1, i])
            v_s = 0.5*(v[j, i] + v[j-1, i])
            u_e = 0.5*(u[j-1, i+1] + u[j, i+1])
            u_w = 0.5*(u[j-1, i] + u[j, i])
            convection = -((rho*u_e*v_e - rho*u_w*v_w)/dx)-((rho*v_n**2 - rho*v_s**2)/dy)
            diffusion = mu*((v[j, i+1]-2*v[j, i]+v[j, i-1])/dx**2 + (v[j+1, i]-2*v[j, i]+v[j-1, i])/dy**2)
            # calculating the intermediate velocity
            v[j, i] += (del_t*(convection+diffusion)/rho) 
    return u, v