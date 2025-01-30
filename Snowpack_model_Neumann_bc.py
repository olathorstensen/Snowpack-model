### 1D snowpack temperature simulator ###
# v2.6 Neumann boundary conditions 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
#   - Radiometer and Tinytag data import
#   - Import include: 
#       ¤ Crop to 1.april period
#       ¤ Noice filtering
#       ¤ Interpolation to 30 sec and 1 min intervals

# Comment: 
#   - TODO Change color for temp plot, fix "neuman cond" in Heat_flow function,
#     add increasing density with depth (maybe), make plot code more compact (itr list)
#   - Instenses to looked at marked with "FIX"
#   - Should fix neagtative SW values

# NOTE: Changes made in hour_angle and T1 amp

import numpy as np
import math as mt
import pandas as pd
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import time
start_time = time.time()

############################    Parameters   ############################ 
runtime = 24*3              # Hours (int)
dt = 30                   # Time step [seconds] (Must be a divisor of 3600) For Sc. 4 use 30s or 60s
dx = 0.005                 # Dist. interval [m]
depth = 1                 # Snow depth from surface [m]
b_bc = 0                  # Bottom boundary condition, fixed [degC]
pisp = 4                  # Plot interval spacer [hours] (int)
plot_depth = 0.35         # Depth shown in plots measured from surface [m]
# Spin up
spin_up = 0               # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*7         # Spin-up run time [Hours]
sp_pisp = 24              # Spin-up Plot interval spacer [hours] (int)
# Data and files
load_ic = 1               # Load IC from file [0] No, [1] Yes
ic_to_file = 0            # Writes model IC to file. If spin-up[1] -> IC given by end of spin-up temp. [0] No, [1] Yes
data_to_file = 1          # Write radiation and atm temp data to new file (spin-up excluded) [0] No, [1] Yes

scenario = 4              # Choose scenario [3],[4] Not fully working yet...
window_size = 30          # Rolling window for radiometer data noice reduction
neumann = 1




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
sw_mod = 659.8259         # Theoretical peak shortwave [W/m2]
sw_rr = sw_mod/sw_peak   # Reduction ratio
sw_k = 100               # Solar extinction coefficient (k) 

    # Sensible/latent/Longwave heat flux properties
sigma = 5.67*10**-8      # Stefan-Boltzmann constant
emis = 1                 # Snow emissivity 
C = 0                    # Coefficient 
A = k/dx/6.655           # A-variable
B = sigma * emis         # B-variable
T1_amp = 5               # T1 amplitude [deg C]  #6
T1_avg = -4              # T1 average temp. [deg C]
T1_phase = 6600          # T1 phase shift [s]

T2_amp = 3               # T2 amplitude [deg C]
T2_avg = -31             # T2 average temp. [deg C]
T2_phase = 6600          # T2 phase shift [s] 



############################    Data import   ############################ 
current_dir = os.path.dirname(os.path.abspath(__file__))  

if load_ic > 0:
    #input_file = r'D:\Dokumenter\...\IC_data_output.xlsx' # Use manual file path if current_dir doesnt work
    input_file = os.path.join(current_dir, 'IC_data.xlsx')
    print("File loaded:", input_file)
    
    df1 = pd.read_excel(input_file, header=0)
    row = df1.iloc[:,0].values
    ic = np.array(row, dtype=float)
    
if scenario  == 4:
    
    # Radiometer data
    input_file = os.path.join(current_dir, 'Radiometer_data.xlsx')
    print("File loaded:", input_file)
    df_r = pd.read_excel(input_file, header=11)
    columns = ['Date and time','Net SW (W/m2)', 'Net LW (W/m2)', 'ELWup', 'ELWlow', 'Snow surface temp (°C)']
    df_r = df_r[columns]
    # Coarse crop and noice reduction
    start_filter = pd.to_datetime("2022-03-29" + " 00:00:00") 
    end_filter = pd.to_datetime("2022-04-03" + " 00:00:00")
    df_r = df_r[df_r["Date and time"].between(start_filter, end_filter)]
    df_r = df_r.reset_index(drop=True)
    df_rf = df_r[["Date and time"]].copy() # New dataframe: 'radiometer_filtered'
    for i in columns[1:]:
        df_rf[i] = df_r[i].rolling(window=window_size, center=True).max() \
                          .rolling(window=window_size, center=True).median() \
                          .rolling(window=window_size, center=True).mean()  
    # Final crop
    start_filter = pd.to_datetime("2022-03-30" + " 00:00:00") 
    end_filter = pd.to_datetime("2022-04-02" + " 00:00:00")
    df_rf = df_rf[df_rf["Date and time"].between(start_filter, end_filter)] # Selecting 1.april period
    df_rf = df_rf.reset_index(drop=True)    
 
    if dt == 30:
        freq = "30S" # Interpolates 30 sec interval
        df_rf.set_index("Date and time", inplace=True)
        df_rf = df_rf.resample(freq).interpolate(method="linear")
        df_rf.reset_index(inplace=True) 
    else:
        freq = "1T"  # Interpolates 1 min interval
        
    # Chose input columns and convert to numpy array with formula below    
    rad_data = np.array(df_rf.iloc[:,5].values, dtype=float)
    SW_net = np.array(df_rf.iloc[:,1].values, dtype=float)#Choose a column FIX
    # Lwin = ....
    # SWnet = ...

    # Tinytag temperature data
    input_file = os.path.join(current_dir, 'Tinytag_data.xlsx')
    print("File loaded:", input_file)
    df_t = pd.read_excel(input_file, header=0)
    df_t = df_t[['Date_B', 'Tinytag_B']] # Tinytag B - 150 cm
    df_t.set_index("Date_B", inplace=True)
    df_t = df_t.resample(freq).interpolate(method="linear")
    df_t.reset_index(inplace=True)   

    df_t = df_t[df_t["Date_B"].between(start_filter, end_filter)] # Selecting 1.april period
    atm_data = np.array(df_t.iloc[:,1].values, dtype=float)
 
    print("atm data shape", atm_data.shape)
    print("rad_data", rad_data.shape)

###########################    Functions    ###########################
def Solar_rad(h_angle, iy):
    elevation = mt.degrees(mt.asin(mt.cos(mt.radians(lat))* mt.cos((mt.radians(h_angle)))))
    zenith = 90 - elevation
    if elevation > 0:
        radiation = sw_con * mt.cos(mt.radians(zenith))
    else:
        radiation = 0   
    rad_scaled = radiation/sw_rr
    if spin_up == 0:
        solar_array[0,iy] = h_angle
        solar_array[1,iy] = elevation
        solar_array[2,iy] = zenith
        solar_array[3,iy] = radiation
        solar_array[4,iy] = rad_scaled 
    return rad_scaled


def Solar_extinction(x0,sw):
    # Need depth and scaled solar radiation
    x0 = x0/100 # con. cm -> m
    x1 = x0+dx
    sw_joules = sw*dt*(1 - mt.exp(-sw_k*depth)) # W/m2*s -> Joules
    sw_partition = ((sw_joules*mt.exp(-sw_k*x0)) - (sw_joules*mt.exp(-sw_k*x1))) /(cp*mass) # J-> dT
    if spin_up == 0:
        solar_array[5,iy] = sw_joules
    return sw_partition


def Heat_flux_surface(t, ix, iy, grid, sw_in) :
    iy = iy-1
    T = grid[ix, iy]        # Temp. surface
    T_dx = grid[ix+1, iy]   # Temp. x= +dx
    T1 = -T1_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T1_phase)) + T1_avg # Temp air
    T2 = -T2_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T2_phase)) + T2_avg # Temp atm   
    #heat_flux = T+(C+A*(T1-T) - B*((T+T0)**4 - (T2+T0)**4))*(dx/k)               
    heat_flux = (
    T_dx 
    + (C + A * (T1 - T) 
        - B * ((T + T0)**4 - (T2 + T0)**4) 
        + ((sw_in * mt.exp(-sw_k * 0)) 
          - (sw_in * mt.exp(-sw_k * dx))))
    * (2*dx / k)
    )
    if spin_up == 0:
        # Writes input data to file apart from during spin-up
        T_array[0,iy] = T1
        T_array[1,iy] = T2
        T_array[2,iy] = heat_flux
        T_array[3,iy] = T
        
    return heat_flux


def Heat_flow(ix, iy, grid, neumann, bc_cell):
    iy = iy-1 
    if neumann == 1 and ix == 0:
            t1 = bc_cell[iy]
    else:                                       
        t1 = grid[ix-1, iy]   
    t2 = grid[ix, iy] 
    t3 = grid[ix+1, iy]
    T = t2 + r*((-2 * t2) + t1 + t3) 
    return T


def Latent_heat(T_value):
    if T_value > 0:
        L = T_value             
    elif T_value <= 0:
        L = 0
    return L
    

def Vapor_pressure(ix, iy, T):
    T = T + 273.15
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
        fgr = 1.1297 * (np.log(v) - np.log(5))
    elif v < -5:
        fgr = -1 * (1.1297 * (np.log(np.abs(v)) - np.log(5)))
    else:
        fgr = 0    
    return fgr


def Facet_growth(ix, iy):
    fg = (fgr[ix, iy] + fgr[ix + 1, iy]) / 2 * dt / 10**6 
    return fg


def Diurnal_array_reshape(array, runtime):
    array_dummy = np.delete(array, 0) 
    if runtime > 24:
        for i in np.arange(1, int(runtime/24)):
            array = np.concatenate((array, array_dummy), axis=0)
        # add remainder
        if runtime%24 != 0:
            array_ext = array_dummy[0:(int(runtime%24*h))]
            array = np.concatenate((array, array_ext))

    if runtime < 24:
        array = array[0:int(runtime*h)+1]  # Crops temp forcing
    return array
    
    
#########################    Model domain    #########################

    # Model grid  
nx = int(depth/dx)
ny = int((runtime * 60 * 60) / dt)

temp = np.zeros([nx+1, ny+1], dtype=float)    # Main snowpack grid
temp_dummy = np.zeros([nx+1, ny+1], dtype=float) 
latent = np.zeros([nx+1, ny+1], dtype=float)  # Latent heat grid
vp = np.zeros([nx+1, ny+1], dtype=float)      # Vapor pressure grid
vpg = np.zeros([nx+1, ny+1], dtype=float)     # Vapor pressure gradient grid
fgr = np.zeros([nx+1, ny+1], dtype=float)     # Facet growth rate grid
fg =  np.zeros([nx+1, ny+1], dtype=float)     # Facet growth grid


    # Axis and time
x = np.linspace(0, depth*100, nx+1) # Depth axis [cm]
y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
y_sec = y.astype(float)
y_hours = np.round(np.arange(0, ny+1, 1) *dt/3600)

base_time = datetime.strptime("00:00 01/01/2025", "%H:%M %d/%m/%Y") 
y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M %d/%m/%Y") for seconds in y_sec]


    # Initial condition
# Linear ic
if load_ic == 0:
    ic = np.linspace(-16, 0, nx +1)    
temp[:, 0] = np.round(ic, 5)


    # Boundary condition
# Bottom BC
temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)

# Surface BC    
ghost_cell = np.zeros([ny+1], dtype=float)
hour_angle = np.linspace(-202, 158, int(h * 24) + 1) #FIX -180, 180
hour_angle = Diurnal_array_reshape(hour_angle, runtime) #Extend or crops input to runtime length
solar_array = np.zeros([6,ny+1], dtype=float) # Content: hour angle, elevation, zenith, rad, rad scaled, joules
T_array = np.zeros([4,ny+1], dtype=float)


# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number:', r)
if r > 0.5:
    print('The r_number is to high, should be < 0.5. Try adjust dt or dx')
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', temp.shape)



#####################     Spin up    ##################### 
if spin_up == 1:
    print('Spin-up initiated. Spin-up time', sp_runtime, 'hours')
    sp_ny = int((sp_runtime * 60 * 60) / dt) #Spin-up time steps
    sp_y = np.arange(0, sp_ny+1, 1) * dt #Spin-up y-axis [sec]
    sp_temp = np.zeros([nx+1, sp_ny+1], dtype=float) # Spin-up temp grid
    sp_latent = np.zeros([nx+1, sp_ny+1], dtype=float) # Spin-up latent heat grid
    sp_ghost_cell = np.zeros([sp_ny+1], dtype=float)
    sp_hour_angle = Diurnal_array_reshape(hour_angle, sp_runtime)
    sp_temp[:,0] = ic
    sp_temp[-1,:] = b_bc * np.ones(sp_ny+1, dtype=float)  #Fixed bottom bc

            
    for iy in np.arange(1, sp_ny+1, dtype=int):       
        for ix in np.arange(0, nx, dtype=int):
            sw_in = Solar_rad(sp_hour_angle[iy], iy)
            if ix == 0:
                #Surface temp calc.
                sp_ghost_cell[iy-1] = Heat_flux_surface(sp_y[iy], ix, iy, sp_temp, sw_in)
                sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + sp_latent[ix, iy-1]
                
            else:
                # Temp. for snow pack
                sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + Solar_extinction(x[ix],sw_in) + sp_latent[ix,iy]
            
            sp_latent[ix,iy] = Latent_heat(sp_temp[ix,iy])
            if sp_temp[ix, iy] > 0:
                sp_temp[ix, iy] = 0        

    spin_up = 0     
    ic = sp_temp[:, -1]
    temp[:,0] = ic

    plt.plot(sp_temp[:,-1], x, label= f"Time {sp_y[-1]/3600} h") # Plot every whole hour    
    plt.gca().invert_yaxis()
    plt.title('Spin-up end state used for IC')
    plt.xlabel('Temperature [C] ')
    plt.ylabel('Depth [cm]')
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()  

    xticks = np.arange(0,-20,-2)
    pld = int(plot_depth/dx)
    for p in np.arange(0, sp_ny+1, h*sp_pisp):
        plt.plot(sp_temp[:pld,p], x[:pld], label= f"{sp_y[p]/3600} Hours") # Plot every whole hour
     
    plt.gca().invert_yaxis()
    plt.title('Scenario 3: Temperature spin-up')
    plt.xlabel('Temperature C')
    plt.ylabel('Depth [cm]')
    plt.legend(fontsize=8)
    plt.grid(alpha=0.5)
    plt.xticks(xticks)
    plt.show()
    



#####################     Main Model Loop    ##################### 
# Temperature

#temp[0,:] = rad_data

for iy in np.arange(1, ny+1, dtype=int):       
    for ix in np.arange(0, nx, dtype=int):  #Change 1 to 0 for neumann
        sw_in = Solar_rad(hour_angle[iy], iy)
        if ix == 0:
            #Surface temp. calc.
            ghost_cell[iy-1] = Heat_flux_surface(y[iy], ix, iy, temp, sw_in)
            temp[ix, iy] = Heat_flow(ix, iy, temp, neumann, ghost_cell) + latent[ix, iy-1] 
        else:
            # Snowpack temp. calc.
            temp[ix, iy] = Heat_flow(ix, iy, temp, neumann, ghost_cell) + Solar_extinction(x[ix],sw_in) + latent[ix, iy-1]           
            #temp[ix, iy] = Heat_flow(ix, iy, temp, neumann, ghost_cell) + latent[ix, iy-1]
            
        latent[ix, iy] = Latent_heat(temp[ix, iy])
        if temp[ix, iy] > 0:
            temp[ix, iy] = 0
        

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



##########################     Plots    ########################## 
xticks = np.arange(0,-20,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx)



for p in np.arange(0, ny+1, h*pisp):
    plt.plot(temp[:pld,p], x[:pld], label= f"{y_t[p]}") # Plot every whole hour
 
plt.gca().invert_yaxis()
plt.title('Scenario 3: Temperature')
plt.xlabel('Temperature C')
plt.ylabel('Depth [cm]')
plt.legend(fontsize=8)
plt.grid(alpha=0.5)
plt.xticks(xticks)
plt.show()


for p in np.arange(0, ny+1, h*pisp):
    plt.plot(vp[:pld,p], x[:pld], label= f"Time {y_t[p]}") # Plot every whole hour

plt.gca().invert_yaxis()
plt.title('Vapor pressure')
plt.xlabel('vp [mb] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


# for p in np.arange(0, ny+1, h*pisp):
#     plt.plot(vpg[:pld,p], x[:pld], label= f"{y_t[p]}") # Plot every whole hour
    
# plt.gca().invert_yaxis()
# plt.title('Vapor pressure gradient')
# plt.xlabel('vpg [mb/m] ')
# plt.ylabel('Depth [cm]')
# plt.legend()
# plt.grid(alpha=0.5)
# plt.show()


# for p in np.arange(0, ny+1, h*pisp):
#     plt.plot(fgr[:pld,p], x[:pld], label= f"{y_t[p]}") # Plot every whole hour
    
# plt.gca().invert_yaxis()
# plt.title('Facet growth rate')
# plt.xlabel('Facet growth rate [nm/s] ')
# plt.ylabel('Depth [cm]')
# plt.legend()
# plt.grid(alpha=0.5)
# plt.show()


# plt.plot(net_growth[:pld], x[:pld])
# plt.gca().invert_yaxis()
# plt.title('Scenario 3: Net facet growth')
# plt.xlabel('Net growth [mm] ')
# plt.ylabel('Depth [cm]')
# #plt.legend()
# plt.grid(alpha=0.5)
# plt.show()



############################    Data output   ############################ 
if data_to_file == 1:
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
    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False)
    print("Wrote to file:", output_file)
    
if ic_to_file == 1:
    #Write end-state to file for IC
    output_file2 = os.path.join(current_dir, 'IC_data.xlsx')
    data_titles2 = [f"Spin-up temp after {sp_y[-1]} hours"]
    data2 = {
        data_titles2[0]: temp[:,0],
        }
    df = pd.DataFrame(data2)
    df.to_excel(output_file2, index=False)
    print("Wrote to file:", output_file2)

end_time = time.time()
print(f"Simulation complete. Runtime: {(end_time-start_time):.2f} seconds")



### Plot-Block: Can be runned independently from model routine in Spyder  
#%%
#Surface temp
plt.plot(y,rad_data, label= "Measured") # Plot every whole hour  
plt.plot(y, sp_temp[0,11520:], label= "Spin-up final 3-days") # Plot every whole hour   
plt.title('Surface temp')
plt.xlabel('Temperature [C] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

#%%
# #Solar plot
# plt.plot(y[:2880], solar_array[4,:2880], label= "Calc") # Plot every whole hour    
# plt.plot(y[:2880], SW_net[:2880], label= "Measured 30.03.22") # Plot every whole hour  
# plt.title('Shortwave rad')
# plt.xlabel('Temperature [C] ')
# plt.ylabel('Depth [cm]')
# plt.legend()
# plt.grid(alpha=0.5)
# plt.show()

#Solar plot
plt.plot(y, solar_array[4,:], label= "Calc") # Plot every whole hour    
plt.plot(y, SW_net, label= "Measured 30.03.22") # Plot every whole hour  
plt.title('Shortwave rad')
plt.xlabel('Temperature [C] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

#%%
#Temperature
plt.plot(y, T_array[0,:], label= "T1") # Plot every whole hour   
plt.plot(y, atm_data, label= "Tinytag B 150cm") # Plot every whole hour  
plt.title('SW in')
plt.xlabel('Seconds')
plt.ylabel('Temperature [C]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

# plt.plot(y, atm_data, label= f"Time {sp_y[-1]/3600} h") # Plot every whole hour    
# plt.gca().invert_yaxis()
# plt.title('SW in')
# plt.xlabel('Temperature [C] ')
# plt.ylabel('Depth [cm]')
# plt.legend()
# plt.grid(alpha=0.5)
# plt.show()
