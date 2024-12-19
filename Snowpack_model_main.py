# -*- coding: utf-8 -*-
### 1D snowpack temperature simulator ###
# v1.3
#@author: olath
# Verison update: VP function

import numpy as np
import matplotlib.pyplot as plt



############################    Parameters   ############################ 
runtime = 24         # Hours
dt = 200                 # Time step [seconds] (Must be devisible by 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.02               # Dist. interval [m]

sin_ini = 0              # 0 for liner ic, 1 for sinus test ic
pm = 4                   # USE INTEGERS. Plot multiplier, 1 -> every h, 2 -> every second h


############################    Constants   ############################   
k = 0.2                  # Thermal conductivity snow [W/m K]
rho = 200                # density [kg/m3]
cp = 2090                # Spesific heat ice [J/kg C]
t0 = 273.15          # Kelvin [K]
a = k/(rho*cp)           # Thermal diffusivity [m2/s]
r = a*(dt/(dx*dx))       # Must be < 0.5 for model stability [1]
h = int((3600/dt))       # Number of dt increments between each whole hour [1]




###########################    Functions    ###########################
def Heat_diff(ix, iy):
    iy = iy-1
    t1 = snow[ix-1, iy]
    t2 = snow[ix, iy] 
    t3 = snow[ix+1, iy]
    
    temp = t2 + r*((-2*t2) + t1+t3)
                  
    return temp


def Vapor_pres(ix, iy, temp):
    temp = temp + t0
    vp = -9.09718 * ((t0/temp)-1) -3.56654 * np.log10(t0/temp) 
    vp = vp +0.876793 * (1-(t0/temp)) + np.log10(6.1071)
    vp = 10**vp
    
    return vp


#########################    Model domain    #########################
 
    # Model grid  
nx = int(depth/dx)
ny = int((runtime*60*60)/dt)
x = np.linspace(0, depth*100, nx+1)
y = np.round( np.arange(0, ny+1, 1) * 200/3600 , 2)

snow = np.zeros([nx+1,ny+1], dtype=float)  # Main snowpack grid
vp = np.zeros([nx+1,ny+1], dtype=float)    # Vapor pressure grid

    ## Initial condition
if sin_ini == 1:
    
    #Sinusoidal test condition
    xx = np.linspace(0, 180, nx+1)
    ic = -4*np.sin(np.deg2rad(xx))-2
    snow[:, 0] = np.round(ic,4)
    
    # Fixed bc
    bc = -2 *np.ones(ny, dtype=float)
    snow[0,:] = bc
    snow[-1,:] = bc
    
else:
    # Linear ic
    ic = np.linspace(-12, -2, nx+1)
    #ic = np.ones(nx+1)*(-2)
    snow[:, 0] = np.round(ic,4)
    
    # Fixed bc, diurnal oscillation:
    xx = np.linspace(-90, 270, int(h*24)+1)
    print('xx shape', xx.shape)
    bc = 5*np.sin(np.deg2rad(xx))-7      #Sinusoidal surface bc
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
        
    
    print('bc shape:', bc.shape)
    print('y', y.shape)
    
    plt.plot(y, bc) 
    plt.title('Temperatrure surface forcing')
    plt.xlabel('Hours')
    plt.show()
    
    snow[0,:] = bc
    snow[-1,:] = -2 * np.ones(ny+1, dtype=float)  #Fixed bottom bc
    
    # Could move forcing param to "Parameter section" 
    
        
# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number:', r)
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', snow.shape)




#####################     Main Model Loop    ##################### 

for iy in np.arange(1, ny+1, dtype=int):
    for ix in np.arange(1, nx, dtype=int):
        snow[ix, iy] = Heat_diff(ix, iy)
        
for iy in np.arange(0, ny+1, dtype=int):
    for ix in np.arange(0, nx+1, dtype=int):
        vp[ix, iy] = Vapor_pres(ix, iy, snow[ix, iy])

###################################################################



# PLot of results:
for p in np.arange(0, ny+1, h*pm):
    plt.plot(snow[:,p], x, label= f"Time {y[p]} h") # Plot every whole hour

    
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


