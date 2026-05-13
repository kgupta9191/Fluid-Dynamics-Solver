import numpy as np
import matplotlib.pyplot as plot

lx = 1.0
ly = 1.0

nx = 8
ny = 8

dx = lx/nx
dy = ly/ny

x = np.linspace(0, lx-dx, nx)
y = np.linspace(0, ly-dy, ny)
X, Y = np.meshgrid(x, y)

swirl_strength = 1.0
core_radius = 0.5

# Calculate radial distance from center
r = np.sqrt((X - 0.5)**2 + (Y - 0.5)**2)

# Generate velocity components for the swirl flow
u = -swirl_strength * (Y - 0.5) * np.exp(-r**2 / (2 * core_radius**2))
v = swirl_strength * (X - 0.5) * np.exp(-r**2 / (2 * core_radius**2))

U = np.zeros((ny, nx))
V = np.zeros((ny, nx))
velMag = np.zeros((ny, nx))

# Interpolate my cell centered U and V velocites
for i in range(0, nx):
    for j in range(0, ny):
        U[j, i] = u[j, i]
        V[j, i] = v[j, i]
        velMag[j, i] = np.sqrt(U[j, i]**2 + V[j, i]**2)


##################################

x_l = np.linspace(0.75*dx, lx-(0.25*dx), nx)
y_l = np.linspace(0.75*dy, ly-(0.25*dy), ny)
X_l, Y_l = np.meshgrid(x_l, y_l)


# Calculate radial distance from center
r_l = np.sqrt((X - 0.5)**2 + (Y - 0.5)**2)

# Generate velocity components for the swirl flow
u_l = -swirl_strength * (Y_l - 0.5) * np.exp(-r_l**2 / (2 * core_radius**2))
v_l = swirl_strength * (X_l - 0.5) * np.exp(-r_l**2 / (2 * core_radius**2))

U_l = np.zeros((ny, nx))
V_l = np.zeros((ny, nx))

# Interpolate my cell centered U and V velocites
for i in range(0, nx):
    for j in range(0, ny):
        U_l[j, i] = u_l[j, i]
        V_l[j, i] = v_l[j, i]


######################################


u_new = np.zeros((ny+4, nx+4))
v_new = np.zeros((ny+4, nx+4))
u_new[2:ny+2, 2:nx+2] = u
v_new [2:ny+2, 2:nx+2] = v 
u_lag = np.zeros((ny+4, nx+4))

#Boundary condition
u_new[0, 0:nx+3] = u_new[3, 0:nx+3]
u_new[1, 0:nx+3] = u_new[2, 0:nx+3]
u_new[ny+3, 0:nx+3] = u_new[ny, 0:nx+3]
u_new[ny+2, 0:nx+3] = u_new[ny+1, 0:nx+3]
u_new[0:ny+3, 0] = u_new[0:ny+3, 3]
u_new[0:ny+3, 1] = u_new[0:ny+3, 2]
u_new[0:ny+3, nx+3] = u_new[0:ny+3, nx]
u_new[0:ny+3, nx+2] = u_new[0:ny+3, nx+1]
v_new[0:ny+3, 0] = v_new[0:ny+3, 3]
v_new[0:ny+3, 1] = v_new[0:ny+3, 2]
v_new[0:ny+3, nx+3] = v_new[0:ny+3, nx]
v_new[0:ny+3, nx+2] = v_new[0:ny+3, nx+1]
v_new[0, 0:nx+3] = v_new[3, 0:nx+3]
v_new[1, 0:nx+3] = v_new[2, 0:nx+3]
v_new[ny+3, 0:nx+3] = v_new[ny, 0:nx+3]
v_new[ny+2, 0:nx+3] = v_new[ny+1, 0:nx+3]

u_rx = np.zeros((ny+4, nx+4))
u_ry = np.zeros((ny+4, nx+4))
v_rx = np.zeros((ny+4, nx+4))
v_ry = np.zeros((ny+4, nx+4))

def phi(r):
     if r <= -2 and r >= 2: 
           return 0
     if r >= -2 and r <= -1: 
           return (5+2*r-np.sqrt(-7-12*r-4*r**2))/8
     if r >= -1 and r <= 0: 
           return (3+2*r+np.sqrt(1-4*r-4*r**2))/8
     if r >= 0 and r <= 1: 
           return (3-2*r+np.sqrt(1+4*r-4*r**2))/8
     if r >= 1 and r <= 2: 
           return (5-2*r-np.sqrt(-7+12*r-4*r**2))/8
     return 0
     

u_x_eul = np.linspace(0-(2*dx), lx+dx, nx+4)
u_y_eul = np.linspace(-1.5*dy, ly+(1.5*dy), ny+4)
v_y_eul = np.linspace(0-(2*dy), ly+dy, ny+4)
v_x_eul = np.linspace(-1.5*dx, lx+(1.5*dx), nx+4)
u_X_eul, u_Y_eul = np.meshgrid(u_x_eul, u_y_eul)
v_X_eul, v_Y_eul = np.meshgrid(v_x_eul, v_y_eul)
X_lag, Y_lag = np.meshgrid(u_x_eul, v_y_eul)
X_lag = (0.75*dx)+X_lag
Y_lag = (0.75*dx)+Y_lag

u_int = np.zeros((ny+4, nx+4))
v_int = np.zeros((ny+4, nx+4))

for i in range(2, ny+2):
            for j in range(2, nx+2): 
                 for f in range (ny+4):
                      for h in range (nx+4):
                             u_rx[h, f] = (X_lag[j,i] - u_X_eul[h,f])/dx 
                             u_ry[h, f] = (Y_lag[j,i] - u_Y_eul[h,f])/dy 
                             u_int[h, f] = u_new[h, f] * phi(u_rx[h, f]) * phi(u_ry[h, f])
                 u_lag[j,i] = np.sum(u_int) 
                 
v_lag = np.zeros((ny+4, nx+4))
for i in range(2, ny+2):
            for j in range(2, nx+2): 
                 for f in range (ny+4):
                      for h in range (nx+4):
                             v_rx[h, f] = (X_lag[j,i] - v_X_eul[h,f])/dx 
                             v_ry[h, f] = (Y_lag[j,i] - v_Y_eul[h,f])/dy 
                             v_int[h, f] = v_new[h, f] * phi(v_rx[h, f]) * phi(v_ry[h, f])
                 v_lag[j,i] = np.sum(v_int) 
                 
U_lag = np.zeros((ny, nx))
V_lag = np.zeros((ny, nx))
velMag_lag = np.zeros((ny, nx))
U_diff = np.zeros((ny, nx))
V_diff = np.zeros((ny, nx))

# Interpolate my cell centered U and V velocites
for i in range(0, nx):
    for j in range(0, ny):
        U_lag[j, i] = u_lag[j+2, i+2] 
        V_lag[j, i] = v_lag[j+2, i+2]
        velMag_lag[j, i] = np.sqrt(U_lag[j, i]**2 + V_lag[j, i]**2)

# Interpolate my cell centered U and V velocites
for i in range(0, nx):
    for j in range(0, ny):
        U_diff[j, i] = U[j, i] - U_lag[j, i]
        V_diff[j, i] = V[j, i] - V_lag[j, i]



# Plot stream particles
fig1, ax1 = plot.subplots(figsize=(6, 6))
strm = ax1.streamplot(x, y, U, V, color=U, linewidth=2, cmap=plot.cm.autumn)
cbar = plot.colorbar(strm.lines)
cbar.ax.set_ylabel('Velocity magnitude')
ax1.set_title('Stream Particles')
ax1.set_xlabel('X-coordinate (m)')
ax1.set_ylabel('Y-coordinate (m)')

# Plot streamlines
#fig3, ax3 = plot.subplots()
#startX, startY = np.meshgrid(np.arange(0, lx, 0.2), np.arange(0, ly, 0.3))
#strm = ax3.streamplot(x, y, U, V, start_points=np.array([startX, startY]).T, linewidth=1.5, arrowsize=1.5)
#ax3.set_title('Velocity Streamline')
#ax3.set_xlabel('X-coordinate (m)')
#ax3.set_ylabel('Y-coordinate (m)')

# Plot velocity vectors
fig2, ax2 = plot.subplots()
q = ax2.quiver(x_l, y_l, U_l, V_l, units='width', color='k', scale=2)
ax2.set_title('Eulerian velocity')
ax2.set_xlabel('X-coordinate (m)')
ax2.set_ylabel('Y-coordinate (m)')

# Plot velocity vectors
fig4, ax4 = plot.subplots()
q = ax4.quiver(x, y, U_lag, V_lag, units='width', color='k', scale=2)
ax4.set_title('Lagrangian velocity')
ax4.set_xlabel('X-coordinate (m)')
ax4.set_ylabel('Y-coordinate (m)')

# Plot velocity vectors
fig5, ax5 = plot.subplots()
q = ax5.quiver(x, y, U_diff, V_diff, units='width', color='k', scale=2)
ax5.set_title('Difference velocity')
ax5.set_xlabel('X-coordinate (m)')
ax5.set_ylabel('Y-coordinate (m)')

plot.show()

