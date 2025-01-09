### 1D snowpack temperature simulator ###
# v2.1 Neumann boundary conditions 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
#   - Can write solar and T1/T2 temp data to file
#   - Can load Simons IC, but pathway must be specified manually to file

# Comment: 
#   - Max possible runtime: 24h
#   - Need to add new vp function


import numpy as np
import math as mt
import pandas
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta


############################    Parameters   ############################ 
runtime = 24             # Hours
dt = 30                  # Time step [seconds] (Must be a divisor of 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.005               # Dist. interval [m]
pm = 2                   # USE INTEGERS. Plot hour interval
b_bc = 0                 # Bottom boundary condition, fixed [degC]
spin_up = 0              # [0] No spin-up, [1] Run spin-up
sp_runtime = 24          # Spin-up run time [Hours]
plot_depth = 0.35        # Depth shown in plots measured from surface [m]
data_to_file = 1         # Write radiation and atm temp data to new file [0] No, [1] Yes
load_ic_file = 1         # Load IC from file [0] No, [1]Steady, [2]+20%, [3]+30%, [4]+70, [5]-20%




############################    Constants   ############################ 

    # Snow properties  
k = 0.2                  # Thermal conductivity snow [W/m K]
rho = 277                # density [kg/m3]
cp = 2090                # Specific heat capacity of ice [J/kg C]
T0 = 273.15              # Ice-point temperature [K]
a = k/(rho*cp)           # Thermal diffusivity [m2/s]
r = a*(dt/(dx*dx))       # Must be < 0.5 for model stability [1]
mass = rho*dx            # Mass of one snow cell with area 1 m2 [Kg]
h = int((3600/dt))       # Number of dt increments between each whole hour [1]

    # Solar properties
lat = 61                 # Latitude [deg]
sw_con = 1361            # Solar constant [W/m2]
sw_peak = 110            # Observed shortwave peak [W/m2]
sw_mod = 659.8259        # Theoretical peak shortwave [W/m2]
sw_rr = sw_mod/sw_peak   # Reduction ratio
sw_k = 100               # Solar extinction coefficient (k)

    # Sensible/latent/Longwave heat flux properties
sigma = 5.67*10**-8      # Stefan-Boltzmann constant
emis = 1                 # Snow emissivity 
C = 0                    # Coefficient 
A = k/dx/6.655           # A-variable
B = sigma * emis         # B-variable
T1_amp = 6               # T1 amplitude [deg C]
T1_avg = -4              # T1 average temp. [deg C]
T1_phase = 6600          # T1 phase shift [s]

T2_amp = 3               # T2 amplitude [deg C]
T2_avg = -31             # T2 average temp. [deg C]
T2_phase = 6600          # T2 phase shift [s]



############################    Data import   ############################ 
    
if load_ic_file > 0:
    input_file = r'D:\Dokumenter\Akademisk arbeid\Near surface Vapour pressure - paper\Snowpack model\IC_input.xlsx'
    #current_dir = os.path.dirname(os.path.abspath(__file__)) #Needs bug fix
    #input_file = os.path.join(current_dir, 'IC_input.xlsx')  #Needs bug fix
    df1 = pandas.read_excel(input_file, header=None)
    row = df1.iloc[load_ic_file-1, 1:].values
    ic = np.array(row, dtype=float)


###########################    Functions    ###########################
def Solar_rad(h_angle, iy):
    elevation = mt.degrees(mt.asin(mt.cos(mt.radians(lat))* mt.cos((mt.radians(h_angle)))))
    solar_array[0,iy] = elevation
    
    zenith = 90 - elevation
    if elevation > 0:
        radiation = sw_con * mt.cos(mt.radians(zenith))
    else:
        radiation = 0   
    rad_scaled = radiation/sw_rr
    solar_array[0,iy] = h_angle
    solar_array[1,iy] = elevation
    solar_array[2,iy] = zenith
    solar_array[3,iy] = radiation
    solar_array[4,iy] = rad_scaled
    return rad_scaled


def Solar_extinction(x0,x1,sw):
    # Need depth and scaled solar radiation
    x0 = x0/100 # con. cm -> m
    x1 = x1/100 # con. cm -> m
    sw_joules = sw*dt*(1-mt.exp(-sw_k*depth)) # Radiation convertet to Joules
    solar_array[5,iy] = sw_joules
    A_sw = sw_joules * sw_k                   # A-term in the SW extinction formula
    sw_partition = ((A_sw/-sw_k) * (mt.exp(-sw_k*x1) - mt.exp(-sw_k*x0))) /(cp*mass)
    return sw_partition


def Sensible_longwave_heat(t, srf_T, iy):
    iy = iy-1
    T1 = -T1_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T1_phase)) + T1_avg # Temp air
    T2 = -T2_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T2_phase)) + T2_avg # Temp atm
    T = srf_T    
    heat_flux = T+(C+A*(T1-T) - B*((T+T0)**4 - (T2+T0)**4))/(k/dx)
    T_array[0,iy] = T1
    T_array[1,iy] = T2
    T_array[2,iy] = heat_flux
    T_array[3,iy] = T
    return heat_flux


def Heat_flow(ix, iy, grid, neuman):
    if neuman == 1:
        if ix == 0:
            iy = iy-1                                        
            t1 = ghost_cell[iy]
            t2 = grid[ix, iy] 
            t3 = grid[ix+1, iy]
            T = t2 + r*((-2 * t2) + t1 + t3)       
        else:
            iy = iy-1                                        
            t1 = grid[ix-1, iy]
            t2 = grid[ix, iy] 
            t3 = grid[ix+1, iy]
            T = t2 + r*((-2 * t2) + t1 + t3)          
    else:
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
x = np.linspace(0, depth*100, nx+1) # Depth axis [cm]
y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
y_sec = y.astype(float)
y_hours = np.round(np.arange(0, ny+1, 1) *dt/3600)

base_time = datetime.strptime("00:00", "%H:%M")
y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M") for seconds in y_sec]


    # Initial condition
# Linear ic
if load_ic_file == 0:
    ic = np.linspace(-16, 0, nx +1)
    
temp[:, 0] = np.round(ic, 5)


    # Boundary condition
# Bottom BC
temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)

# Surface BC    
ghost_cell = np.zeros([ny+1], dtype=float)
hour_angle = np.linspace(-180, 180, int(h * 24) + 1)
solar_array = np.zeros([6,ny+1], dtype=float) # Content: hour angle, elevation, zenith, rad, rad scaled, joules
T_array = np.zeros([4,ny+1], dtype=float)

   
# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number:', r)
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', temp.shape)



#####################     Spin up    ##################### 

xx = np.linspace(-90, 270, int(h * 24) + 1)
bc = 9*np.sin(np.deg2rad(xx))-10      #Sinusoidal surface bc
bc_dummy = np.delete(bc, 0)
sp_bc = bc

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
            sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 0)

          
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
# Temperature
for iy in np.arange(1, ny+1, dtype=int):
    ghost_cell[iy-1] = Sensible_longwave_heat(y[iy], temp[0,iy-1], iy)
        
    for ix in np.arange(0, nx, dtype=int):    
        sw_in = Solar_rad(hour_angle[iy], iy)
        temp[ix, iy] = Heat_flow(ix, iy, temp, 1) + Solar_extinction(x[ix+0],x[ix+1],sw_in) 
        

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
xticks = np.arange(0,-20,-2)
pd = int(plot_depth/dx)

for p in np.arange(0, ny+1, h*pm):
    plt.plot(temp[:pd,p], x[:pd], label= f"{y_t[p]}") # Plot every whole hour
 
plt.gca().invert_yaxis()
plt.title('Scenario 3: Temperature')
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
plt.title('Scenario 3: Net facet growth')
plt.xlabel('Net growth [mm] ')
plt.ylabel('Depth [cm]')
#plt.legend()
plt.grid(alpha=0.5)
plt.show()



############################    Data output   ############################ 
if data_to_file == 1:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = os.path.join(current_dir, 'SC3_data_output.xlsx')
    data_titles = ["Time", "Surf_T", "T1", "T2","Hour Angle", "Elevation", 
                   "Zenith", "Radiation", "Rad_scaled", "Rad_en_joules", "A-term sw"]
    data = {
        data_titles[0]: y_t,
        data_titles[1]: T_array[3,:],
        data_titles[2]: T_array[0,:],
        data_titles[3]: T_array[1,:],
        data_titles[4]: solar_array[0,:],
        data_titles[5]: solar_array[1,:],
        data_titles[6]: solar_array[2,:],
        data_titles[7]: solar_array[3,:],
        data_titles[8]: solar_array[4,:],
        data_titles[9]: solar_array[5,:],
        }
    df = pandas.DataFrame(data)
    df.to_excel(output_file, index=False)
