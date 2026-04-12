import numpy as np
import matplotlib.pyplot as plot
import constant 

# Domain Setup
# Length of domain
lx = constant.Geometry_in_x_direction
ly = constant.Geometry_in_y_direction
# No of cells
nx = constant.Cells_in_x_direction
ny = constant.Cells_in_y_direction
# Cell size
dx = lx/nx
dy = ly/ny
# Cell cavity
nxc = constant.Cavity_cells_in_x_direction
nyc = constant.Cavity_cells_in_y_direction

x = np.linspace(0, lx, nx)
y = np.linspace(0, lx, ny)
X, Y = np.meshgrid(x, y)

# Fluid properties
rho = constant.Density # Density 
mu = constant.Viscosity # Viscosity
nu = mu/rho

# Initialize vaiables
u = np.zeros((ny+2,nx+2)) # Face centered
v = np.zeros((ny+2,nx+2)) # Face centered
p = np.zeros((ny+2,nx+2)) # Cell centered
div = np.zeros((ny+2,nx+2)) # Cell centered

u_bc = np.zeros((ny+2,nx+2))  # Initial velocities
v_bc = np.zeros((ny+2,nx+2))

intialVel = constant.Initial_Velocity # Velocity of lid
del_t = min ((0.25*dx*dy)/nu, 4*nu/intialVel**2) 
twfin = 100 * del_t

# Boundary conditions
# Top wall
u_t = intialVel   # Slip
v_t = 0.0         # Penetration
# Bottom wall
u_b = 0.0
v_b = 0.0
# Left wall
u_l = 0.0         # Penetration
v_l = 0.0         # Slip
# Right wall
u_r = 0.0
v_r = 0.0
# Cavity wall
u_c = 0.0
v_c = 0.0

myFile = open('pressureSolverConvergenceData.txt', 'w')
myFile.write('{:>8s} {:>8s} {:>8s}\n'.format('Timestep', 'Residual', '#-of-iterations'))