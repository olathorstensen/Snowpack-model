# -*- coding: utf-8 -*-
### 1D snowpack temperature simulator ###
# v1.6
#@author: Ola Thorstensen and Thor Parmentier
# Version update: Språkvask


import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots



############################    Parameters   ############################ 
runtime = 24         # Hours
dt = 120                 # Time step [seconds] (Must be a divisor of 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.01               # Dist. interval [m]
pm = 6                   # USE INTEGERS. Plot multiplier, 1 -> every h, 2 -> every second h


############################    Constants   ############################   
k = 0.1                  # Thermal conductivity snow [W/m K]
rho = 200                # density [kg/m3]
cp = 2090                # Specific heat capacity of ice [J/kg C]
T0 = 273.16              # Ice-point temperature [K]
a = k/(rho*cp)           # Thermal diffusivity [m2/s]
r = a*(dt/(dx*dx))       # Must be < 0.5 for model stability [1]
h = int((3600/dt))       # Number of dt increments between each whole hour [1]




###########################    Functions    ###########################
def Heat_flow(ix, iy):
    iy = iy-1                                        
    t1 = temp[ix-1, iy]
    t2 = temp[ix, iy] 
    t3 = temp[ix+1, iy]
    
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
x = np.linspace(0, depth*100, nx+1)
y = np.round( np.arange(0, ny+1, 1) * 200/3600 , 2)                             # Should this be dt/3600?

temp = np.zeros([nx+1, ny+1], dtype=float)  # Main snowpack grid
vp = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure grid
vpg = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure gradient grid
fgr = np.zeros([nx+1, ny+1], dtype=float)    # Facet growth rate grid
fg =  np.zeros([nx+1, ny+1], dtype=float)    # Facet growth grid

    ## Initial condition
# Linear ic
ic = np.linspace(-10, 0, nx+1)
temp[:, 0] = np.round(ic, 4)

# Fixed bc, diurnal oscillation:
xx = np.linspace(-90, 270, int(h * 24) + 1)
bc = 9*np.sin(np.deg2rad(xx))-10      #Sinusoidal surface bc
bc_dummy = np.delete(bc, 0)

# Extends or crops 24h temperature swings to runtime length
if runtime > 24:
    for i in np.arange(1, int(runtime/24)):
        bc = np.concatenate((bc, bc_dummy,), axis=0)
    # add remainder
    if runtime%24 != 0:
        bc_dummy = bc_dummy[0:(int(runtime%24*h))]
        bc = np.concatenate((bc, bc_dummy))

if runtime < 24:
    bc = bc[0:int(runtime*h)+1]  # Crops temp forcing

plt.plot(y, bc) 
plt.title('Temperature surface forcing')
plt.xlabel('Hours')
plt.show()

temp[0,:] = bc
temp[-1,:] = 0 * np.ones(ny+1, dtype=float)  #Fixed bottom bc

# Could move forcing param to "Parameter section" 
    
        
# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number:', r)
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', temp.shape)




#####################     Main Model Loop    ##################### 

for iy in np.arange(1, ny+1, dtype=int):
    for ix in np.arange(1, nx, dtype=int):
        temp[ix, iy] = Heat_flow(ix, iy)
        
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
for p in np.arange(0, ny+1, h*pm):
    plt.plot(temp[:,p], x, label= f"Time {y[p]} h") # Plot every whole hour
 
plt.gca().invert_yaxis()
plt.title('Temperature')
plt.xlabel('Temperature C')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(vp[:,p], x, label= f"Time {y[p]} h") # Plot every whole hour

plt.gca().invert_yaxis()
plt.title('Vapor pressure')
plt.xlabel('vp [mb] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(vpg[:,p], x, label= f"Time {y[p]} h") # Plot every whole hour
    
plt.gca().invert_yaxis()
plt.title('Vapor pressure gradient')
plt.xlabel('vpg [mb/m] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


for p in np.arange(0, ny+1, h*pm):
    plt.plot(fgr[:,p], x, label= f"Time {y[p]} h") # Plot every whole hour
    
plt.gca().invert_yaxis()
plt.title('Facet growth rate')
plt.xlabel('Facet growth rate [nm/s] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

plt.plot(net_growth, x)
plt.gca().invert_yaxis()
plt.title('Net facet growth')
plt.xlabel('Net growth [mm] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()



