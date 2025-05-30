### 1D Snowpack Temperature SimOlaThor ###
# v3.9 Scenario 3 and 4 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
# - Small work on plotting sc3 and 4
# - Should have spinup more than 20 days


# Comment: 
#   - Instenses to looked at marked with "FIX"

#Storagee of latent heat

#FIX Solar curve peak at 13:36
# Scale theoretical solar input

import numpy as np
import math as mt
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from datetime import datetime, timedelta
import time

start_time = time.time()

############################    Parameters   ############################ 
scenario = 3                # Must be 3 or 4. Read paper for explanation, dummy.                  
dt = 30                    # Time step [seconds] (Must be a divisor of 3600) For Sc. 4 use 30s or 60s
dx = 0.005                 # Dist. interval [m]
depth = 1                  # Snow depth from surface [m]
b_bc = 0                   # Bottom boundary condition, fixed [°C]
pisp = 2                   # Plot interval spacer [hours] (int)
plot_depth = 0.35          # Depth shown in plots measured from surface [m]

spin_up = 1               # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*21        # Spin-up run time [Hours]
sp_pisp = 24              # Spin-up Plot interval spacer [hours] (int)

load_ic = 0               # Load IC from file [0] No, [1] Yes
ic_to_file = 0            # Writes model IC to file. If spin-up[1] -> IC given by end of spin-up temp. [0] No, [1] Yes
data_to_file = 0          # Write radiation and atm temp data to new file (spin-up excluded) [0] No, [1] Yes
ic_scaling = 1            # How warm or cold the IC should be. Warmer IC  [0.7], colder IC [1.3]


cold_T = 1                # For SC3, use [0] for org T temp, use [1] for 700m incresed elevation
window_size = 30          # Rolling window for radiometer data noise reduction
sw_k = 50                # Solar extinction coefficient (k) 

ng_title = 'K = 50 $m^{-1}$'           # Title for DTVPGE dataframe
ng_switch = 1



############################    Parameters multi runs   ############################ 
if scenario == 3:
    runs = 1                 # Number of model runs (spinup + main model) must be
    ng_title_list = [ng_title,
                     'K = 30 $m^{-1}$',
                     'K = 100 $m^{-1}$',
                     'Warmer IC',
                     'Colder IC']           # Title for dtvpge dataframe
    sw_k_list = [sw_k, 30, 100, 50, 50]                                   # Solar extinction coefficient (k)
    ic_scaling_list = [ic_scaling, 1, 1, 0.7, 1.3]                  # Scales IC Warmer IC [0.7], Colder IC [1.3]
    runtime = 24*2                                              # Hours (int)
    skew_sw = 0                                               # For sc 3 or maybe only 4? I think only 4
    
if scenario == 4:
    runs = 1                     # Number of model runs (spinup + main model) must be
    runtime = 24*3            # Hours (int)
    skew_sw = 1              # For sc 3 or maybe only 4? I think only 4
    
    
        
############################    Constants   ############################ 

    # Snow properties  
k = 0.1439275              # Thermal conductivity snow [W/m K]    0.1439
rho = 245                # density [kg/m3] #277
cp = 2090                # Specific heat capacity of ice [J/kg °C]
T0 = 273.15              # Ice-point temperature [K]
a = k/(rho*cp)           # Thermal diffusivity [m2/s]
r = a*(dt/(dx*dx))       # Must be < 0.5 for model stability [1]
mass = rho*dx            # Mass of one snow cell with area 1 m2 [kg]
h = int((3600/dt))       # Number of dt increments between each whole hour [1]

    # Solar properties
lat = 61                 # Latitude [°]
sw_con = 1361            # Solar constant [W/m2]
sw_peak = 110            # Observed shortwave peak [W/m2] 
sw_mod = 659.8259        # Theoretical peak shortwave [W/m2]
sw_rr = sw_mod/sw_peak   # Reduction ratio


    # Sensible/latent/Longwave heat flux properties
sigma = 5.67*10**-8      # Stefan-Boltzmann constant
emis = 1                 # Snow emissivity 
C = 0                    # Coefficient 
A = (k/dx)/4.14         # A-variable #3.8  #6.655 org SC3 4.14
B = sigma * emis         # B-variable
T1_amp = 6               # T1 amplitude [°C]  #6
T1_avg = -4              # T1 average temp. [°C]
T1_phase = 6600          # T1 phase shift [s]

T2_amp = 3               # T2 amplitude [°C]
T2_avg = -31             # T2 average temp. [°C]
T2_phase = 6600          # T2 phase shift [s] 


# roling window: 1.30 , 2. 100, 3. 100
############################    Data import   ############################ 
current_dir = os.path.dirname(os.path.abspath(__file__))  

# #sIMON TEMPS
# input_file = os.path.join(current_dir, 'surface_temp_simon.xlxs')
# print("File loaded:", input_file)
# df_simon = pd.read_excel(input_file, header=0)

# SW scaling data - Works only for dt=30s
input_file = os.path.join(current_dir, 'CR674_scaling.csv')
print("File loaded:", input_file)
df_s = pd.read_csv(input_file, header=0)
SW_scaling = np.array(df_s.iloc[:,0].values, dtype=float)
        

if load_ic > 0:
    #input_file = r'D:\Dokumenter\...\IC_data_output.xlsx' # Use manual file path if current_dir doesnt work
    input_file = os.path.join(current_dir, 'IC_data.xlsx')
    print("File loaded:", input_file)
    
    df1 = pd.read_excel(input_file, header=0)
    row = df1.iloc[:,0].values
    ic = np.array(row, dtype=float)

if scenario == 3:
    input_file = os.path.join(current_dir, 'DTVPGE_data.xlsx')
    print("File loaded:", input_file)
    df_ng = pd.read_excel(input_file, header=0)    
    
if scenario  == 4:
    # Radiometer data
    input_file = os.path.join(current_dir, 'Radiometer_data.xlsx')
    print("File loaded:", input_file)
    df_r = pd.read_excel(input_file, header=11)
    columns = ['Date and time','Net SW', 'Net LW', 'ELWup', 'ELWlow', 'Snow surface temp']
    df_r = df_r[columns]
    # Coarse crop and noise reduction
    start_filter = pd.to_datetime("2022.03.29" + " 00:00:00") 
    end_filter = pd.to_datetime("2022.04.03" + " 00:00:00")
    df_r = df_r[df_r["Date and time"].between(start_filter, end_filter)]
    df_r = df_r.reset_index(drop=True)
    df_rf = df_r[["Date and time"]].copy() # New dataframe: 'radiometer_filtered'
    for i in columns[1:]:
        df_rf[i] = df_r[i].rolling(window=window_size, center=True).max() \
                          .rolling(window=window_size, center=True).median() \
                          .rolling(window=window_size, center=True).mean()  
    # Final crop
    start_filter = pd.to_datetime("2022.03.30" + " 00:00:00") 
    end_filter = pd.to_datetime("2022.04.02" + " 00:00:00")
    df_rf = df_rf[df_rf["Date and time"].between(start_filter, end_filter)] # Selecting 1.april period
    df_rf = df_rf.reset_index(drop=True)    
    # org df
    df_r = df_r[df_r["Date and time"].between(start_filter, end_filter)] # Selecting 1.april period
    df_r = df_r.reset_index(drop=True)    
 
 
    if dt == 30:
        freq = "30S" # Interpolates 30 sec interval
        df_rf.set_index("Date and time", inplace=True)
        df_rf = df_rf.resample(freq).interpolate(method="linear")
        df_rf.reset_index(inplace=True) 
        #org df
        df_r.set_index("Date and time", inplace=True)
        df_r = df_r.resample(freq).interpolate(method="linear")
        df_r.reset_index(inplace=True) 
    else:
        freq = "1T"  # Interpolates 1 min interval  
    # Convert from pandas to numpy
    rad_data = df_rf['Snow surface temp'].values.astype(float)
    # SW_net = np.array(df_rf.iloc[:,1].values, dtype=float)
    # SW_net = np.where(SW_net<0, 0, SW_net) #Removes negative values
    # LW_net = np.array(df_rf.iloc[:,2].values, dtype=float)
    SW_net = np.array(df_r.iloc[:,1].values, dtype=float)
    SW_net = np.where(SW_net<0, 0, SW_net) #Removes negative values
    LW_net = np.array(df_r.iloc[:,2].values, dtype=float)
  

    # Tinytag temperature data
    input_file = os.path.join(current_dir, 'Tinytag_data.xlsx')
    print("File loaded:", input_file)
    df_t = pd.read_excel(input_file, header=0)
    for i in ['A','B']:
        tt_date = 'Date_'+str(i)
        tt_temp = 'Tinytag_'+str(i)
        df_tf = df_t[[tt_date, tt_temp]] 
        df_tf.set_index(tt_date, inplace=True)
        df_tf = df_tf.resample(freq).interpolate(method="linear")
        df_tf.reset_index(inplace=True)   
        df_tf = df_tf[df_tf[tt_date].between(start_filter, end_filter)] # Selecting 1.april period
        if i == 'A':  
            tinytag_A = np.array(df_tf.iloc[:,1].values, dtype=float)
        elif i == 'B':
            tinytag_B = np.array(df_tf.iloc[:,1].values, dtype=float)
            
    
    SW_scaled = SW_net * SW_scaling


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
    if cold_T == 0:
        T1 = -T1_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T1_phase)) + T1_avg # Temp air org
        T2 = -T2_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T2_phase)) + T2_avg # Temp atm org 
    else:
        lr = (0.0065)*700 # laps rate -6.5 deg/1000m, up 'x' m
        T1 = -T1_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T1_phase)) + (T1_avg-lr) # Temp air colder, 
        T2 = -T2_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T2_phase)) + (T2_avg-lr) #0.36)) # Temp atm colder, 36% reduction in laps rate
    
    Q_sensible = C + A*(T1-T)
    Q_LWout = -B * ((T + T0)**4 - (T2 + T0)**4) 
    
    heat_flux = (
    T_dx 
    + (Q_sensible
     + Q_LWout 
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
        T_array[4,iy] = Q_sensible
        T_array[5,iy] = Q_LWout
        
    return heat_flux


def Heat_flow(ix, iy, grid, bc_type, bc_cell):
    iy = iy-1 
    if bc_type == 1 and ix == 0:
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
    fg = (fgr[ix, iy] + fgr[ix, iy+1]) / 2 * dt / 10**6 
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
 


#########################    Model loop    #########################   
for run in range (1, runs + 1): 
    if runs > 1:
        ng_title = ng_title_list[run - 1]
        sw_k = sw_k_list[run - 1]
        ic_scaling = ic_scaling_list[run - 1]
    
    
    #########################    Model domain    #########################
            
        # Model grid  
    nx = int(depth/dx)
    ny = int((runtime * 60 * 60) / dt)
    
    temp = np.zeros([nx+1, ny+1], dtype=float)    # Main snowpack grid
    temp_dummy = np.zeros([nx+1, ny+1], dtype=float) 
    latent = np.zeros([nx+1, ny+1], dtype=float)  # Latent heat grid
    vp = np.zeros([nx+1, ny+1], dtype=float)      # Vapor pressure grid
    vpg = np.zeros([nx, ny+1], dtype=float)     # Vapor pressure gradient grid
    fgr = np.zeros([nx, ny+1], dtype=float)     # Facet growth rate grid
    fg =  np.zeros([nx, ny], dtype=float)     # Facet growth grid
    dtvpge = np.zeros_like(vpg, dtype=float)
    
    
        # Axis and time
    x = np.linspace(0, depth*100, nx+1) # Depth axis [cm]
    x_stag = x[:-1]+dx*100/2            # Staggered x axis for vpg, fgr, fg, ng, and dtvpge
    y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
    y_sec = y.astype(float)
    y_hours = np.round(np.arange(0, ny+1, 1) *dt)/3600
    
    base_time = datetime.strptime("00:00 30/03/2022", "%H:%M %d/%m/%Y") 
    y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M %d/%m/%Y") for seconds in y_sec]
    
    
    # Initial condition
    # Linear IC
    if load_ic == 0:
        ic = np.linspace(-16, 0, nx +1)   
    temp[:, 0] = ic
    
    # Boundary conditions
    # Bottom BC
    temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)
    
    # Surface BC    
    ghost_cell = np.zeros([ny+1], dtype=float)
    hour_angle = np.linspace(-204, 156, int(h * 24) + 1)
    hour_angle = Diurnal_array_reshape(hour_angle, runtime) # Extend or crops input to runtime length
    solar_array = np.zeros([6,ny+1], dtype=float) # Content: hour angle, elevation, zenith, rad, rad scaled, joules
    T_array = np.zeros([6,ny+1], dtype=float)
    
    
    # Prints of shapes and numbers (all are numbers and none are shapes)
    print('r_number (should be less than 0.5):', r)
    if r > 0.5:
        print('The r_number is too high, should be < 0.5. Try adjusting dt or dx')


    
    #####################     Spin up    ##################### 
    spin_up_has_occurred = 0
    
    if spin_up == 1:
        print('Spin-up initiated. Spin-up time', sp_runtime, 'hours')
        sp_ny = int((sp_runtime * 60 * 60) / dt) # Spin-up time steps
        sp_y = np.arange(0, sp_ny+1, 1) * dt # Spin-up y-axis [sec]
        sp_temp = np.zeros([nx+1, sp_ny+1], dtype=float) # Spin-up temp grid
        sp_latent = np.zeros([nx+1, sp_ny+1], dtype=float) # Spin-up latent heat grid
        sp_ghost_cell = np.zeros([sp_ny+1], dtype=float)
        sp_hour_angle = Diurnal_array_reshape(hour_angle, sp_runtime)
        sp_temp[:,0] = ic
        sp_temp[-1,:] = b_bc * np.ones(sp_ny+1, dtype=float)  # Fixed bottom bc
    
        if skew_sw == 1:
            sp_SW_scaling = Diurnal_array_reshape(SW_scaling[0:2881], sp_runtime)
           
        for iy in np.arange(1, sp_ny+1, dtype=int):       
            for ix in np.arange(0, nx, dtype=int):
                sw_in = Solar_rad(sp_hour_angle[iy], iy)
                if skew_sw == 1:
                    sw_in = sw_in * sp_SW_scaling[iy]
                if ix == 0:
                    #Surface temp calc
                    sp_ghost_cell[iy-1] = Heat_flux_surface(sp_y[iy], ix, iy, sp_temp, sw_in)
                    sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + sp_latent[ix, iy-1]
                    
                else:
                    # Temp for snowpack
                    sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + Solar_extinction(x[ix],sw_in) + sp_latent[ix,iy-1]
                
                sp_latent[ix,iy] = Latent_heat(sp_temp[ix,iy])
                if sp_temp[ix, iy] > 0:
                    sp_temp[ix, iy] = 0    
                 
        spin_up_has_occurred = 1    
        ic = sp_temp[:, -1]
        temp[:,0] = ic
        
        # Indexes the last 24h of the spin up surface temp
        index = int((86400/dt)*((sp_runtime/24)-1))
        sp_temp_mean = np.mean(sp_temp[0,index:-1])
        print('Diurnal mean surface temperature', sp_temp_mean) 
        print('Spin-up complete')
        
        
    
    #####################     Main Model Calculations    ##################### 
    # Temperature
    
    if scenario == 4:
        bc_type = 0           # [0] Dirichlet (fixed)
        linscale = rad_data[0]/ic[0]
        #linscale = tinytag_A[0]/ic[20]
        temp[0,:] = rad_data      # Surface BC
        temp[:,0] =  temp[:,0] * linscale # Scale IC
        #temp[:,0] = temp[:,0] * ic_scaling
        
    elif scenario == 3:
        bc_type = 1         # [1] Neumann (ghost cell)
        #temp = sp_temp[:,]
        temp[:,0] = temp[:,0] * ic_scaling
    
    
    for iy in np.arange(1, ny+1, dtype=int):       
        for ix in np.arange(1 if bc_type == 0 else 0, nx, dtype=int):
            
            if scenario == 3:
                sw_in = Solar_rad(hour_angle[iy], iy)
                if skew_sw == 1:
                    sw_in = sw_in * SW_scaling[iy]
                if ix == 0:
                    #Surface temp. calc.
                    ghost_cell[iy-1] = Heat_flux_surface(y[iy], ix, iy, temp, sw_in)
                    temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + latent[ix, iy-1] 
                else:
                    # Snowpack temp. calc.
                    temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + Solar_extinction(x[ix],sw_in) + latent[ix, iy-1]  
                    
            elif scenario == 4:
                Solar_rad(hour_angle[iy], iy)
                if skew_sw == 1:
                    temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + Solar_extinction(x[ix],SW_scaled[iy]) + latent[ix, iy-1]
                else:
                    temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + Solar_extinction(x[ix],SW_net[iy]) + latent[ix, iy-1]
                  
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
            
    for iy in np.arange(0, ny, dtype=int):
        for ix in np.arange(0, nx, dtype=int):
            fg[ix, iy] = Facet_growth(ix, iy)
    
    if scenario == 4:
        net_growth = np.sum(fg[:,h*24:h*48], axis = 1)
        SC4_net_growth = np.zeros([3,int(depth/dx)])
        SC4_net_growth[0,:] = np.sum(fg[:,h*0:h*24], axis = 1)
        SC4_net_growth[1,:] = np.sum(fg[:,h*24:h*48], axis = 1)
        SC4_net_growth[2,:] = np.sum(fg[:,h*48:h*72], axis = 1)
        
    else:
        #net_growth = np.sum(fg, axis = 1)
        net_growth = np.sum(fg[:,0:h*24], axis = 1)
    
    # DTVPGE calculations
    dtvpge = np.where(vpg < -5, vpg + 5, np.where(vpg > 5, vpg - 5, 0))
   
    if scenario == 4:
        dtvpge_mean = np.mean(dtvpge[:,h*24:h*48], axis = 1)
        SC4_dtvpge_mean = np.zeros([3,int(depth/dx)])
        SC4_dtvpge_mean[0,:] = np.mean(dtvpge[:,h*0:h*24], axis = 1)
        SC4_dtvpge_mean[1,:] = np.mean(dtvpge[:,h*24:h*48], axis = 1)
        SC4_dtvpge_mean[2,:] = np.mean(dtvpge[:,h*48:h*72], axis = 1)
    else:
        dtvpge_mean = np.mean(dtvpge, axis = 1)
    
    if run == 1:
        dtvpge_abs = np.abs(dtvpge)
        dtvpge_abs_mean = np.mean(dtvpge_abs, axis = 1)
        
    
    
    ############################    Data Output   ############################ 
    # Sensible heat calculation
    if scenario == 3:
        if run == 1:
            SC3_sensible_heat = T_array[4,:]
            SC3_LW_out = T_array[5,:]
            SC3_srf_temp = temp[0,:]
            SC3_10cm_temp = temp[20,:]
            SC3_net_growth = net_growth
            SC3_dtvpge_mean = dtvpge_mean
        if run == 4:
            SC3_dtvpge_mean_Warmer_IC = dtvpge_mean
    elif scenario == 4:
        sw_srf = SW_scaled*(1 - mt.exp(-sw_k*dx))
        sensible_heat = (((temp[0,:]-temp[1,:])/dx)*k) - sw_srf - LW_net
        
        LW_out = LW_net # FIX this should be changed
    
    if data_to_file == 1:
        output_file = os.path.join(current_dir, 'SC3_data_output.xlsx')
        data_titles = ["Time", "Surf_T", "T1", "T2", "Sensible heat", 
                       "LW out", "Hour Angle", "Elevation", 
                       "Zenith", "Radiation", "Rad_scaled", "Rad_en_joules", "A-term sw"]
        data = {
            data_titles[0]: y_t,
            data_titles[1]: T_array[3,:],
            data_titles[2]: T_array[0,:],
            data_titles[3]: T_array[1,:],
            data_titles[4]: T_array[4,:],
            data_titles[5]: T_array[5,:],
            data_titles[6]: solar_array[0,:],
            data_titles[7]: solar_array[1,:],
            data_titles[8]: solar_array[2,:],
            data_titles[9]: solar_array[3,:],
            data_titles[10]: solar_array[4,:],
            data_titles[11]: solar_array[5,:],
            }
        df = pd.DataFrame(data)
        df.to_excel(output_file, index=False)
        print("Wrote to file:", output_file)
        
    if ic_to_file == 1 and spin_up_has_occurred == 1:
        #Write end-state to file for IC
        output_file2 = os.path.join(current_dir, 'IC_data.xlsx')
        data_titles2 = [f"Spin-up temp after {sp_y[-1]} hours"]
        data2 = {
            data_titles2[0]: sp_temp[:,-1],
            }
        df = pd.DataFrame(data2)
        df.to_excel(output_file2, index=False)
        print("Wrote to file:", output_file2)
     
    
    # output_file2 = os.path.join(current_dir, 'ng_hub_SC3.xlsx')
    # ng_hub.to_excel(output_file2, index=False)
    
    # DTVPGE to df
    if ng_switch == 1:
        if run == 1 and scenario == 3:
            ng_data1 = {ng_title: dtvpge_mean,}
            ng_df1 = pd.DataFrame(ng_data1)
            ng_hub = pd.concat([df_ng, ng_df1], axis=1)
        elif run == 1 and scenario == 4:
            ng_data1 = {ng_title: dtvpge_mean,}
            ng_df1 = pd.DataFrame(ng_data1)
            ng_data0 = {'SC 3': SC3_dtvpge_mean}
            df_3 = pd.DataFrame(ng_data0)
            ng_hub = pd.concat([df_3, ng_df1], axis=1)
            
        else:    
            ng_data2 = {ng_title: dtvpge_mean}
            ng_df2 = pd.DataFrame(ng_data2)
            ng_hub = pd.concat([ng_hub, ng_df2], axis=1)
        print('ng_hub', ng_hub.columns)
    
    if run == 1:
        temp_plot_SC3 = temp                         # When doing multi runs, plot temps from first run


####################################################################################
end_time = time.time()
print(f"Simulation complete. Runtime: {(end_time-start_time):.2f} seconds")


#%%
##########################     Plots    ########################## 


### Plot-Block: Can be run independently from model routine in Spoider  
#%%
yticks = np.arange(120,-140,-20)
if scenario == 3:
    #Solar plot
    plt.plot(y,  solar_array[4,:], label= "Calc SW net", linestyle='--' ) # Plot every whole hour 
    #plt.plot(y, SW_net, label= "Measured SW net") # Plot every whole hour 
    plt.plot(y, LW_out, label= "Calc LW out") # Plot every whole hour 
    plt.plot(y, sensible_heat, label= "Calc sensible heat") # Plot every whole hour
    plt.title('Shortwave rad')
    plt.xlabel('sec')
    plt.ylabel('W/m2')
    plt.yticks(yticks)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()
else:
    #Solar plot
    plt.plot(y, solar_array[4,:], label= "Calc SW net", linestyle='--', color='blue' ) # Plot every whole hour    
    plt.plot(y, SW_net, label= "Measured SW net", linestyle='-', color='blue' ) # Plot every whole hour 
    plt.plot(y, SW_scaled, label= "Scaled SW net") # Plot every whole hour
    plt.plot(y, LW_net, label= "Net LW", linestyle='-', color='green' ) # Plot every whole hour 
    plt.plot(y, SC3_LW_out, label= "Steady state LW out", linestyle='--', color='green' ) # Plot every whole hour 
    plt.plot(y, sensible_heat, label= "Calc sensible heat", linestyle='-', color='orange') 
    plt.plot(y, SC3_sensible_heat, label= "Steady state sensible h.", linestyle='--', color='orange') 
    plt.title('Fluxes') 
    plt.xlabel('sec')
    plt.ylabel('W/m2')
    plt.yticks(yticks)
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()
 

#%%
#Temperature

# fig, ax = plt.subplots(figsize=(10, 6))  # Create figure and axis
# lw = 3.5
# ax.plot(y_t, temp[20, :], label="Snow 10cm simulated", color='C0', lw=lw) 
# ax.plot(y_t, temp[22, :], label="Snow 10cm simulated", color='orange')  # Uncomment to plot this line
# ax.plot(y_t, tinytag_A, label="Tinytag at 10cm", color='C3', lw=lw)  
# ax.plot(y_t, rad_data, label="Surface temperature smoothed", linestyle='-', color='C2', lw=lw)  
# ax.set_title(f'Tinytag vs snow temp. SW k={sw_k}, with latent heat and scaled SW')
# ax.set_xlabel('Time')
# ax.set_ylabel('Temperature [°C]')
# ax.legend()
# ax.grid(alpha=0.5)




# Assuming y_t is a pandas Series or DataFrame column with datetime strings
y_tt = pd.to_datetime(y_t, dayfirst=True)  # Convert y_t to pandas datetime if it's not already

# Convert to pydatetime


import matplotlib.dates as mdates
# Create figure and axis
fig, ax = plt.subplots(figsize=(10, 6))

# Line width
lw = 3

# Plot time series data
ax.plot(y_tt, rad_data, label="Surface temperature smoothed", linestyle='-', color='gold', lw=lw) 
ax.plot(y_tt, SC3_srf_temp , label="Surface SC3 warm", linestyle='--', color='C2', lw=lw) 
ax.plot(y_tt, tinytag_A, label="Tinytag at 10cm", color='C3', lw=lw)  
ax.plot(y_tt, temp[20,:], label="Snow 10cm simulated", color='C0', lw=lw) 


 



# import matplotlib.ticker as ticker
# ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=12))  # Adjust nbins as needed



# Set the date format for the x-axis labels
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M %d-%m-%y'))

# Set the frequency of the x-ticks
ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))  # Show every 1 hour, adjust as necessary
ax.set_xlim([y_tt[0].replace(hour=0, minute=0, second=0), y_tt[-1]])

#ax.set_title(f'Tinytag vs snow temp. SW k={sw_k}, with latent heat and scaled SW')
ax.set_xlabel('Time')
ax.set_ylabel('Temperature [°C]')
plt.xticks(rotation=45)
ax.legend()
ax.grid(alpha=0.5)
fig.tight_layout()

# Show plot




#%%
# Slide show mania

xticks = np.arange(0,-22,-2) 
pld = int(plot_depth/dx)


# Temp slider
fig1, ax1 = plt.subplots(figsize=(9, 6))
plt.subplots_adjust(bottom=0.2)  # Space for slider
line1, = ax1.plot(temp[:pld, 0], x[:pld], label=f'Time: {y_t[0]}')

ax1.set_title("Temperature Profile")
ax1.set_xlabel("Temperature [°C]")
ax1.set_ylabel("Depth [cm]")
ax1.set_xticks(xticks)
ax1.invert_yaxis()
ax1.legend()
ax1.grid(alpha=0.5)

# Slider
slider_ax1 = plt.axes([0.2, 0.05, 0.6, 0.03])  # [left, bottom, width, height]
time_slider = Slider(slider_ax1, "Time Step", 0, ny - 1, valinit=0, valstep=1)

def update_time(val):
    time_idx = int(time_slider.val)
    line1.set_xdata(temp[:pld, time_idx])
    ax1.legend([f'Time: {y_t[time_idx]}'], loc='upper right')
    # fig1.canvas.draw_idle()  # Redraw the figure. Might not actually need this.

time_slider.on_changed(update_time)
plt.show()


# Growth rate slider
fig2, ax2 = plt.subplots(figsize=(9, 6))
plt.subplots_adjust(bottom=0.3)  # Space for slider
line2, = ax2.plot(y_hours, fgr[0, :], label = f'Depth: {x[0]}')

ax2.set_title("DTVPGE")
ax2.set_xlabel("Time [h]")
ax2.set_ylabel("DTVPGE [Pa/cm]")
ax2.set_yticks(np.arange(-150, 200, 50))
ax2.legend()
ax2.grid(alpha=0.5)

# Slider
slider_ax2 = plt.axes([0.2, 0.17, 0.6, 0.03])  # [left, bottom, width, height]
depth_slider = Slider(slider_ax2, "Depth", 0, nx - 1, valinit=0, valstep=1)

# Buttons
button_ax2_up = plt.axes([0.4, 0.1, 0.1, 0.05])  # [left, bottom, width, height]
button_ax2_down = plt.axes([0.5, 0.1, 0.1, 0.05])

button_up = Button(button_ax2_up, "↑ Higher")
button_down = Button(button_ax2_down, "↓ Deeper")

def update_depth(val):
    depth_idx = int(depth_slider.val)
    line2.set_ydata(dtvpge[depth_idx, :])  # Update y-data (FGR at selected depth)
    ax2.legend([f'Depth: {x[depth_idx]:.1f} cm'], loc='upper right')  
    fig2.canvas.draw_idle()  # Redraw the figure. Might not need this

def go_deeper(event):
    current_val = depth_slider.val
    if current_val < nx - 1:
        depth_slider.set_val(current_val + 1)

def go_higher(event):
    current_val = depth_slider.val
    if current_val > 0:
        depth_slider.set_val(current_val - 1)

depth_slider.on_changed(update_depth)
button_down.on_clicked(go_deeper)
button_up.on_clicked(go_higher)

plt.show()


#%%

# PAPER FIGURE SC3
pld = int(plot_depth/dx)+1
plt.rcParams.update({'font.size': 12})
fig, ax = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.1})

# Temperature plot
cmap = plt.get_cmap("turbo_r")
colors = [cmap(i / 12) for i in range(12)]
colors = colors[-5:] + colors [:-5]


for i, p in enumerate(np.arange(0, h*24, h*pisp)):
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax[0].plot(temp_plot_SC3[:pld, p], x[:pld], label=time_only, color=colors[i], lw=2.7)

ax[0].set_xlabel('Temperature [°C]')
ax[0].yaxis.tick_right() 
ax[0].yaxis.set_label_position("right")  
ax[0].set_ylabel('Depth [cm]', rotation=270, labelpad=25)
ax[0].invert_yaxis()
ax[0].legend()
ax[0].grid(alpha=0.5)
ax[0].set_xticks(np.arange(0,-22,-2))
ax[0].xaxis.set_label_position('top')
ax[0].xaxis.tick_top()
ax[0].text(0.968, 0.05, "a",  
           transform=ax[0].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

# Mean DTVPGE near surface
linestyle = ['dotted', '-', '--', '-.', '-', '-' ]
colors = ['C3', 'black', 'darkgrey', 'darkgrey','C1', 'C9' ]
linewidth = [4,4.2,2.7 ,2.7 ,2.7 ,2.7 ,2.7]
for i, column in enumerate(ng_hub.columns):
    ax[1].plot(ng_hub[column][:pld], x_stag[:pld], label=column, linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
  
ax[1].legend()
ax[1].set_xlabel("Mean DTVPGE [Pa/cm]")
ax[1].set_yticklabels([])
ax[1].invert_yaxis()
ax[1].grid(alpha=0.5)
ax[1].xaxis.set_label_position('top')
ax[1].xaxis.tick_top()
ax[1].text(0.93, 0.05, "b",  
           transform=ax[1].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

#%%

# Alternate figures for Scenario 3. Fig 5A6. c plot with zoomed in dtvpge including mean of absolute
plt.rcParams.update({'font.size': 22}) # FIX THIS HAS TO CHANGE
pld = int(plot_depth/dx)+1
fig, ax = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={'width_ratios': [2, 1, 1], 'wspace': 0.1})

# Temperature plot
cmap = plt.get_cmap("turbo_r")
colors = [cmap(i / 12) for i in range(12)]
colors = colors[-5:] + colors [:-5]


for i, p in enumerate(np.arange(0, h*24, h*pisp)):
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax[0].plot(temp_plot_SC3[:pld, p], x[:pld], label=time_only, color=colors[i], lw=2.7)

ax[0].set_xlabel('Temperature [°C]')
ax[0].yaxis.tick_left()
ax[0].set_ylabel('Depth [cm]')
ax[0].yaxis.set_label_position("left")
# ax[0].set_ylabel('Depth [cm]', rotation=270, labelpad=25)
ax[0].invert_yaxis()
ax[0].legend()
ax[0].grid(alpha=0.5)
ax[0].set_xticks(np.arange(0,-22,-2))
ax[0].xaxis.set_label_position('top')
ax[0].xaxis.tick_top()
ax[0].text(0.958, 0.05, "a",  
           transform=ax[0].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

# Mean DTVPGE near surface
linestyle = ['dotted', '-', '--', '-.', '-', '-' ]
colors = ['C3', 'black', 'darkgrey', 'darkgrey','C1', 'C9' ]
linewidth = [4, 4.2, 2.7, 2.7, 2.7, 2.7, 2.7]
for i, column in enumerate(ng_hub.columns):
    ax[1].plot(ng_hub[column][:pld], x_stag[:pld], label=column, linestyle=linestyle[i], color=colors[i], lw=linewidth[i])

ax[1].plot(-1 * dtvpge_abs_mean[:pld], x_stag[:pld], label = 'Mean of absolute', linestyle = 'dotted', color = 'black', lw = 4)
ax[1].legend(loc = 'lower left')
ax[1].set_xlabel("Mean DTVPGE [Pa/cm]")
# ax[1].set_xticks(np.arange(-0.2, 0.06, 0.05))
ax[1].invert_yaxis()
ax[1].grid(alpha=0.5)
ax[1].xaxis.set_label_position('top')
ax[1].xaxis.tick_top()
ax[1].text(0.915, 0.05, "b",  
           transform=ax[1].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

for i, column in enumerate(ng_hub.columns):
    ax[2].plot(ng_hub[column][:pld], x_stag[:pld], label=column, linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
    
ax[2].plot(-1 * dtvpge_abs_mean[:pld], x_stag[:pld], label = 'DTVPGE abs mean', linestyle = 'dotted', color = 'black', lw = 4)
# ax[2].plot(ng_hub['K = 30 $m^{-1}$'][:pld], x_stag[:pld], lw = 2.7, color = 'darkgrey', linestyle = '-.', label = 'K = 30 $m^{-1}$')
# ax[2].plot(ng_hub['K = 50 $m^{-1}$'][0:pld], x_stag[:pld], color = 'black', lw = 2.7, label = 'K = 50 $m^{-1}$')
# ax[2].plot(dtvpge_mean_1[:pld], x_stag[:pld], lw = 2.7, label = '1')
# ax[2].plot(dtvpge_mean_3[:pld], x_stag[:pld], lw = 2.7, label = '3')
# ax[2].plot(dtvpge_mean_14[:pld], x_stag[:pld], lw = 2.7, label = '14')
# # ax[2].plot(dtvpge_mean_70[:pld], x_stag[:pld], lw = 2.7, label = '70')
# ax[2].plot(dtvpge_mean_140[:pld], x_stag[:pld], lw = 2.7, label = '140')
ax[2].set_xlabel("Mean DTVPGE [Pa/cm]")
ax[2].set_xticks(np.arange(-4, 6, 2))
ax[2].set_xlim(-5, 5)
# ax[2].legend(loc='lower left')
ax[2].yaxis.tick_right() 
ax[2].set_ylabel('Depth [cm]')
ax[2].yaxis.set_label_position("right")
ax[2].invert_yaxis()
ax[2].grid(alpha=0.5)
ax[2].xaxis.set_label_position('top')
ax[2].xaxis.tick_top()
ax[2].text(0.915, 0.05, "c",  
           transform=ax[2].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
#%%

# PAPER FIGURE SC4
# To plot SC4, first run SC3 with cold_T = 0, then run SC4 


pld = int(plot_depth/dx)+1
plt.rcParams.update({'font.size': 22})
fig, ax = plt.subplots(1, 3, figsize=(18, 6), gridspec_kw={'width_ratios': [2, 1, 1], 'wspace': 0.1})

#Temperature plot
cmap = plt.get_cmap("turbo_r")
colors = [cmap(i / 12) for i in range(12)]
colors = colors[-5:] + colors [:-5]

for i, p in enumerate(np.arange(h*24, h*48, h*pisp)):
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax[0].plot(temp[:pld, p], x[:pld], label=time_only, color=colors[i], lw=2.7)

ax[0].set_xlabel('Temperature [°C]')
ax[0].yaxis.tick_left()
ax[0].yaxis.set_label_position("left")
ax[0].set_ylabel('Depth [cm]')
# ax[0].set_ylabel('Depth [cm]', rotation=270, labelpad=25)
ax[0].invert_yaxis()
ax[0].legend()
ax[0].grid(alpha=0.5)
ax[0].set_xticks(np.arange(0,-18,-2))
ax[0].xaxis.set_label_position('top')
ax[0].xaxis.tick_top()
ax[0].text(0.96, 0.05, "a", 
           transform=ax[0].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

# Mean DTVPGE near surface

label = ['30.03.22', '31.03.22', '01.04.22']
linestyle = ['-', '-', '-', '--', '--' ]
colors = ['forestgreen', 'limegreen', 'lawngreen', 'black', 'C1' ]
linewidth = [6.5 ,6.5 ,6.5 ,2.4 ,2.4]
for i in np.arange(0,5,1):
    if i == 3:
        ax[1].plot(ng_hub['SC 3'][:pld], x_stag[:pld], label='Scenario 3*', linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
    elif i == 4:
        ax[1].plot(SC3_dtvpge_mean_Warmer_IC[:pld], x_stag[:pld], label='SC3* Warmer IC', linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
    else:
        ax[1].plot(SC4_dtvpge_mean[i,:pld], x_stag[:pld], label=label[i], linestyle=linestyle[i], color=colors[i], lw=linewidth[i])

ax[1].legend()
ax[1].set_xlabel("Mean DTVPGE [Pa/cm]")
ax[1].invert_yaxis()
ax[1].grid(alpha=0.5)
ax[1].xaxis.set_label_position('top')
ax[1].xaxis.tick_top()
ax[1].text(0.91, 0.05, "b",  
           transform=ax[1].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))


for i in np.arange(0,5,1):
    if i == 3:
        ax[2].plot(ng_hub['SC 3'][:pld], x_stag[:pld], label='Scenario 3*$', linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
    elif i == 4:
        ax[2].plot(SC3_dtvpge_mean_Warmer_IC[:pld], x_stag[:pld], label='SC3* Warmer IC', linestyle=linestyle[i], color=colors[i], lw=linewidth[i])
    else:
        ax[2].plot(SC4_dtvpge_mean[i,:pld], x_stag[:pld], label=label[i], linestyle=linestyle[i], color=colors[i], lw=linewidth[i])

ax[2].set_xlabel("Mean DTVPGE [Pa/cm]")
ax[2].set_xlim(-5, 2.5)
# ax[2].legend(loc='lower left')
ax[2].yaxis.tick_right() 
ax[2].set_ylabel('Depth [cm]')
ax[2].yaxis.set_label_position("right")
ax[2].invert_yaxis()
ax[2].grid(alpha=0.5)
ax[2].xaxis.set_label_position('top')
ax[2].xaxis.tick_top()
ax[2].text(0.915, 0.05, "c",  
           transform=ax[2].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
#%%
#Spin up temp profiles
fig, ax = plt.subplots(figsize = (9, 6))

#if spin_up_has_occurred == 1:

# for p in np.arange(0, sp_ny+1, h*sp_pisp):
#     ax.plot(sp_temp[:, p], x, label= f"{sp_y[p]/3600} Hours")
# for p in np.arange(0, ny+1, h*sp_pisp):
#     ax.plot(temp[:, p], x, label= f"{sp_y[p]/3600} Hours")
#ax.plot(temp[:, 0], x, label= f"{y[p]/3600} Hours", color='C1')

ax.plot(sp_y[:2880*3], sp_temp[0, 2880*25:2880*28], label= 'SPIN UP', color='C1')
ax.plot(sp_y[:2880*3], temp[0,:-1], label= 'srf temp', color='C2')
ax.set_title('Temperature profiles spin-up')
ax.set_xlabel('Temperature °C')
ax.set_ylabel('Depth [cm]')
ax.invert_yaxis()
ax.legend(fontsize=8)
ax.grid(alpha=0.5)
ax.set_xticks(xticks)

#%%
# FIGURE 7


plt.figure(figsize=(9, 6))
plt.rcParams.update({'font.size': 14})
plt.plot(y_hours, dtvpge[0,:], label='0 - 5 mm, Surface', lw=2)
plt.plot(y_hours, dtvpge[3,:], label='15 - 20 mm, shows how fast dtvpge decreases', lw=2)
# plt.plot(y_hours, dtvpge[4,:], label='20 - 25 mm, first depth positive dtvpge appears', lw=2)
plt.plot(y_hours, dtvpge[5,:], label='25 - 30 mm', lw=2)
plt.plot(y_hours, dtvpge[9,:], label='45 - 50 mm, first crossing of 0', lw=2)
plt.plot(y_hours, dtvpge[14,:], label='70 - 75 mm, positive peak', lw=2)
plt.plot(y_hours, dtvpge[24,:], label='120 - 125 mm, second crossing of 0', lw=2)
plt.plot(y_hours, dtvpge[36,:], label='180 - 185 mm, lower negative peak', lw=2)
# plt.plot(y_hours, dtvpge[59,:], label='295 - 300 mm, below here only 0', lw=2)
plt.title('')
plt.xlabel('Time [hours')
plt.ylabel('DTVPGE')
plt.xlabel('Time [hours]')
plt.legend(loc="lower left")
plt.grid(True, alpha=0.5)
plt.xticks([0,12,24,36,48])
plt.tight_layout()
plt.show()

#%%
plt.figure(figsize = (9,6))
plt.imshow(dtvpge, #vmin = -25, vmax = 25,
           cmap = 'RdYlGn_r', aspect='auto')
plt.colorbar(label = 'Pa/cm')
plt.title('DTVPGE', size = '16')
plt.xlabel('Time')
plt.ylabel('Depth')

#%%
plt.figure(figsize = (9,6))
plt.imshow(temp, cmap = 'RdBu_r', aspect='auto')
plt.colorbar(label = '°C')
plt.title('Temp', size = '16')
plt.xlabel('Time')
plt.ylabel('Depth')
