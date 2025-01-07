### 1D snowpack temperature simulator ###
# v1.8 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
# - Inplemented parameter for plot depths
# - Time stamp for plots intervals 
# - Fixed error in spin-up module regarding unsymmetrical  model runtime and spin-up runtime
# Comment: Model simulates Scenario 2, but Net facet graph is off at the surface.


import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import plotly.graph_objects as go #Whats this?
import plotly.express as px       #Whats this?
from plotly.subplots import make_subplots #Whats this?



############################    Parameters   ############################ 
runtime = 24             # Hours
dt = 120                 # Time step [seconds] (Must be a divisor of 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.01                # Dist. interval [m]
pm = 2                   # USE INTEGERS. Give hourly plot rate
b_bc = 0                 # Bottom boundary condition, fixed [C]

spin_up = 1            # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*7        # Spin-up run time [Hours]
plot_depth = 0.35           # Depth shown in plots [m]


############################    Constants   ############################   
k = 0.1                  # Thermal conductivity snow [W/m K]
rho = 200                # density [kg/m3]
cp = 2090                # Specific heat capacity of ice [J/kg C]
T0 = 273.16              # Ice-point temperature [K]
a = k/(rho*cp)           # Thermal diffusivity [m2/s]
r = a*(dt/(dx*dx))       # Must be < 0.5 for model stability [1]
h = int((3600/dt))       # Number of dt increments between each whole hour [1]




###########################    Functions    ###########################
def Heat_flow(ix, iy, grid):
    iy = iy-1                                        
    t1 = grid[ix-1, iy]
    t2 = grid[ix, iy] 
    t3 = grid[ix+1, iy]
    T = t2 + r*((-2 * t2) + t1 + t3)     
      
    return T


def Vapor_pressure(ix, iy, T):
    T = T + T0
    vp = -9.09718 * ((T0 / T) - 1) -3.56654 * np.log10(T0 / T) 
    vp = vp +0.876793 * (1-(T/T0)) + np.log10(6.1071)
    vp = 10**vp
    
    return vp

def Vapor_pressure_gradient(ix, iy):
    vpg = (vp[ix, iy] - vp[ix + 1, iy]) / dx 
    
    return vpg

def Facet_growth_rate(ix, iy):
    v = vpg[ix, iy]
    if v > 5:
        fgr = 1.0987 * np.log(v) - 1.7452 
    elif v < -5:
        fgr = -1 * (1.0987 * np.log(np.abs(v)) - 1.7452)
    else:
        fgr = 0 
        
    return fgr


def Facet_growth(ix, iy):
    fg = (fgr[ix, iy] + fgr[ix + 1, iy]) / 2 * dt / 10**6
    
    return fg
    
    
#########################    Model domain    #########################

    # Model grid  
nx = int(depth/dx)
ny = int((runtime * 60 * 60) / dt)

temp = np.zeros([nx+1, ny+1], dtype=float)  # Main snowpack grid
vp = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure grid
vpg = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure gradient grid
fgr = np.zeros([nx+1, ny+1], dtype=float)    # Facet growth rate grid
fg =  np.zeros([nx+1, ny+1], dtype=float)    # Facet growth grid

    # Axis and time
x = np.linspace(0, depth*100, nx+1) # Depth axis
y = np.round( np.arange(0, ny+1, 1) * dt/3600 , 2) # Time axis
base_time = datetime.strptime("00:00", "%H:%M")
y_t = [(base_time + timedelta(hours=hour)).strftime("%H:%M") for hour in y]


    ## Initial condition
# Linear ic
ic = np.linspace(-10, 0, nx+1)
temp[:, 0] = np.round(ic, 4)

# Fixed bc, diurnal oscillation:
xx = np.linspace(-90, 270, int(h * 24) + 1)
bc = 9*np.sin(np.deg2rad(xx))-10      #Sinusoidal surface bc
bc_dummy = np.delete(bc, 0)
sp_bc = bc

# Extends or crops 24h temperature swings to runtime length
if runtime > 24:
    for i in np.arange(1, int(runtime/24)):
        bc = np.concatenate((bc, bc_dummy,), axis=0)
    # add remainder
    if runtime%24 != 0:
        bc_ext = bc_dummy[0:(int(runtime%24*h))]
        bc = np.concatenate((bc, bc_ext))

if runtime < 24:
    bc = bc[0:int(runtime*h)+1]  # Crops temp forcing

plt.plot(y, bc) 
plt.title('Temperature surface forcing')
plt.xlabel('Hours')
plt.show()

temp[0,:] = bc
temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)  #Fixed bottom bc

# Could move forcing param to "Parameter section" 
    
        
# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number:', r)
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', temp.shape)



#####################     Spin up    ##################### 

if spin_up == 1:
    print('Spin-up initiated. Runtime', sp_runtime, 'hours')
    
    if sp_runtime > 24:
        for i in np.arange(1, int(sp_runtime/24)):
            sp_bc = np.concatenate((sp_bc, bc_dummy,), axis=0)
        # add remainder
        if sp_runtime%24 != 0:
            sp_bc_ext = bc_dummy[0:(int(sp_runtime%24*h))]
            sp_bc = np.concatenate((sp_bc, sp_bc_ext))

    if sp_runtime < 24:
        sp_bc = bc[0:int(sp_runtime*h)+1]  # Crops temp forcing


    sp_ny = int((sp_runtime * 60 * 60) / dt) #Spin-up time steps
    sp_y = np.round( np.arange(0, sp_ny+1, 1) * dt/3600 , 2) #Spin-up y-axis
    sp_temp = np.zeros([nx+1, sp_ny+1], dtype=float) # Spin-up temp grid
    sp_temp[:,0] = ic
    sp_temp[0,:] = sp_bc
    sp_temp[-1,:] = b_bc * np.ones(sp_ny+1, dtype=float)  #Fixed bottom bc
    print('dim spin', sp_bc.shape, sp_temp.shape)

    for iy in np.arange(1, sp_ny+1, dtype=int):
        for ix in np.arange(1, nx, dtype=int):
            sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp)

          
    ic = sp_temp[:, -1]
    temp[:,0] = ic

    plt.plot(sp_temp[:,-1], x, label= f"Time {sp_y[-1]} h") # Plot every whole hour    
    plt.gca().invert_yaxis()
    plt.title('Spin-up end state used for IC')
    plt.xlabel('Temperature [C] ')
    plt.ylabel('Depth [cm]')
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()      



#####################     Main Model Loop    ##################### 

for iy in np.arange(1, ny+1, dtype=int):
    for ix in np.arange(1, nx, dtype=int):
        temp[ix, iy] = Heat_flow(ix, iy, temp)
        
for iy in np.arange(0, ny+1, dtype=int):
    for ix in np.arange(0, nx+1, dtype=int):
        vp[ix, iy] = Vapor_pressure(ix, iy, temp[ix, iy])
        
for iy in np.arange(0, ny+1, dtype=int):
    for ix in np.arange(0, nx, dtype=int):
        vpg[ix, iy] = Vapor_pressure_gradient(ix, iy)
        
for iy in np.arange(0, ny+1, dtype=int):
    for ix in np.arange(0, nx, dtype=int):
        fgr[ix, iy] = Facet_growth_rate(ix, iy)
        
for iy in np.arange(0, ny+1, dtype=int):
    for ix in np.arange(0, nx, dtype=int):
        fg[ix, iy] = Facet_growth(ix, iy)

net_growth = np.sum(fg, axis = 1)

###################################################################



# PLot of results:
xticks = np.arange(0,-22,-2)
pd = int(plot_depth/dx)

for p in np.arange(0, ny+1, h*pm):
    plt.plot(temp[:pd,p], x[:pd], label= f"{y_t[p]}") # Plot every whole hour
 
plt.gca().invert_yaxis()
plt.title('Temperature')
plt.xlabel('Temperature C')
plt.ylabel('Depth [cm]')
plt.legend(fontsize=8)
plt.grid(alpha=0.5)
plt.xticks(xticks)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(vp[:pd,p], x[:pd], label= f"Time {y_t[p]}") # Plot every whole hour

plt.gca().invert_yaxis()
plt.title('Vapor pressure')
plt.xlabel('vp [mb] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(vpg[:pd,p], x[:pd], label= f"{y_t[p]}") # Plot every whole hour
    
plt.gca().invert_yaxis()
plt.title('Vapor pressure gradient')
plt.xlabel('vpg [mb/m] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(fgr[:pd,p], x[:pd], label= f"{y_t[p]}") # Plot every whole hour
    
plt.gca().invert_yaxis()
plt.title('Facet growth rate')
plt.xlabel('Facet growth rate [nm/s] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


plt.plot(net_growth[:pd], x[:pd])
plt.gca().invert_yaxis()
plt.title('Net facet growth')
plt.xlabel('Net growth [mm] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()



