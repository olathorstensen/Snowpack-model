### 1D snowpack temperature simulator ###
# v3.2 Scenario 3 and 4 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
#   - T1 and T2 lapse rate scaling
#   - SW scaling (SC4)
#   - Sensbile heat calc (SC3 & SC4)
#   - Use measured SW in SC4





# Comment: 
#   - TODO Change color for temp plot,
#   - Instenses to looked at marked with "FIX"
#   - ic_to file can not be 1 if spin_up is 0
#   - Need to plot the x axis 0.5dx off for vpg, fgr, fg and net_growth
#   - Need staggered grid in y direction for net_growth
#   - Switching scenario 4 data files to csv for increased read speed



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
runtime = 24*3             # Hours (int)
dt = 30                   # Time step [seconds] (Must be a divisor of 3600) For Sc. 4 use 30s or 60s
dx = 0.005                # Dist. interval [m]
depth = 1                 # Snow depth from surface [m]
b_bc = 0                  # Bottom boundary condition, fixed [°C]
pisp = 2                  # Plot interval spacer [hours] (int)
plot_depth = 0.35         # Depth shown in plots measured from surface [m]

spin_up = 1              # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*7        # Spin-up run time [Hours]
sp_pisp = 24              # Spin-up Plot interval spacer [hours] (int)

load_ic = 0               # Load IC from file [0] No, [1] Yes
ic_to_file = 1            # Writes model IC to file. If spin-up[1] -> IC given by end of spin-up temp. [0] No, [1] Yes
data_to_file = 0          # Write radiation and atm temp data to new file (spin-up excluded) [0] No, [1] Yes

bc_type = 0               # [0] Dirichlet (fixed), [1] Neumann (ghost cell)
scenario = 4              # Choose scenario [3],[4] Not fully working yet...
cold_T = 1                # For SC3, use [0] for org T temp, use [1] for 700m incresed elevation
window_size = 30          # Rolling window for radiometer data noice reduction

ng_title = 'IC'           # Title for Net_growth dataframe
run_number = 1            # Must conduct one inital run with [1]
ic_scaling = 1            # Scales IC e.g. 0.7 = -30% scaling

############################    Constants   ############################ 

    # Snow properties  
k = 0.1439              # Thermal conductivity snow [W/m K]    0.1439
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
sw_k = 50                # Solar extinction coefficient (k) 

    # Sensible/latent/Longwave heat flux properties
sigma = 5.67*10**-8      # Stefan-Boltzmann constant
emis = 1                 # Snow emissivity 
C = 0                    # Coefficient 
A = k/dx/4.15             # A-variable #3.8  #6.655 org SC3
B = sigma * emis         # B-variable
T1_amp = 6               # T1 amplitude [°C]  #6
T1_avg = -4              # T1 average temp. [°C]
T1_phase = 6600          # T1 phase shift [s]

T2_amp = 3               # T2 amplitude [°C]
T2_avg = -31             # T2 average temp. [°C]
T2_phase = 6600          # T2 phase shift [s] 



############################    Data import   ############################ 
current_dir = os.path.dirname(os.path.abspath(__file__))  
if run_number < 2:
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
     
        if dt == 30:
            freq = "30S" # Interpolates 30 sec interval
            df_rf.set_index("Date and time", inplace=True)
            df_rf = df_rf.resample(freq).interpolate(method="linear")
            df_rf.reset_index(inplace=True) 
        else:
            freq = "1T"  # Interpolates 1 min interval  
        # Convert from pandas to numpy
        rad_data = np.array(df_rf.iloc[:,5].values, dtype=float)
        SW_net = np.array(df_rf.iloc[:,1].values, dtype=float)
        SW_net = np.where(SW_net<0, 0, SW_net) #Removes negative values
        LW_net = np.array(df_rf.iloc[:,2].values, dtype=float)
  
    
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
                
                
        # SW scaling data - Works only for dt=30s
        input_file = os.path.join(current_dir, 'CR674_scaling.csv')
        print("File loaded:", input_file)
        df_s = pd.read_csv(input_file, header=0)
        SW_scaling = np.array(df_s.iloc[:,0].values, dtype=float)
        
        SW_scaled = SW_net * SW_scaling

        
else:
    print('Variable run_number > 1, Radiometer and Tinytag data not loaded')

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
        T2 = -T2_amp * mt.cos( (2*mt.pi/(24*60*60)) * (t-T2_phase)) + (T2_avg-(lr*0.36)) # Temp atm colder, 36% reduction in laps rate
    
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


    # Axis and time
x = np.linspace(0, depth*100, nx+1) # Depth axis [cm]
x_stag = x[:-1]+dx*100/2            # Staggered x axis for vpg, fgr, fg, and ng
y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
y_sec = y.astype(float)
y_hours = np.round(np.arange(0, ny+1, 1) *dt)/3600

base_time = datetime.strptime("00:00 30/03/2022", "%H:%M %d/%m/%Y") 
y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M %d/%m/%Y") for seconds in y_sec]


    # Initial condition
# Linear ic
if load_ic == 0:
    ic = np.linspace(-16, 0, nx +1)   
temp[:, 0] = ic



    # Boundary conditions
# Bottom BC
temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)

# Surface BC    
ghost_cell = np.zeros([ny+1], dtype=float)
hour_angle = np.linspace(-204, 156, int(h * 24) + 1)
hour_angle = Diurnal_array_reshape(hour_angle, runtime) #Extend or crops input to runtime length
solar_array = np.zeros([6,ny+1], dtype=float) # Content: hour angle, elevation, zenith, rad, rad scaled, joules
T_array = np.zeros([6,ny+1], dtype=float)


# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number (should be less than 0.5):', r)
if r > 0.5:
    print('The r_number is too high, should be < 0.5. Try adjusting dt or dx')
print('a_number:', a)
print('h_number:', h)
print('x', x.shape)
print('y', y.shape)
print('snowpack grid shape', temp.shape)



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

            
    for iy in np.arange(1, sp_ny+1, dtype=int):       
        for ix in np.arange(0, nx, dtype=int):
            sw_in = Solar_rad(sp_hour_angle[iy], iy)
            if ix == 0:
                #Surface temp calc
                sp_ghost_cell[iy-1] = Heat_flux_surface(sp_y[iy], ix, iy, sp_temp, sw_in)
                sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + sp_latent[ix, iy-1]
                
            else:
                # Temp for snowpack
                sp_temp[ix, iy] = Heat_flow(ix, iy, sp_temp, 1, sp_ghost_cell) + Solar_extinction(x[ix],sw_in) + sp_latent[ix,iy]
            
            sp_latent[ix,iy] = Latent_heat(sp_temp[ix,iy])
            if sp_temp[ix, iy] > 0:
                sp_temp[ix, iy] = 0        

    spin_up = 0 
    spin_up_has_occurred = 1    
    ic = sp_temp[:, -1]
    temp[:,0] = ic
    # Indexes the last 24h of the spin up surface temp
    index = int((86400/dt)*((sp_runtime/24)-1))
    sp_temp_mean = np.mean(sp_temp[0,index:-1])
    print('Diurnal mean surface temperature', sp_temp_mean) 
    print('Spin-up complete')
    
    

#####################     Main Model Loop    ##################### 
# Temperature

if scenario == 4:
    #linscale = rad_data[0]/ic[0]
    linscale = tinytag_A[0]/ic[20]
    temp[0,:] = rad_data      # Surface BC
    temp[:,0] = ic * linscale # Scale IC
    temp[:,0] = temp[:,0] * ic_scaling
    
elif scenario == 3:
    bc_type = 1
    temp[:,0] = temp[:,0] * ic_scaling


for iy in np.arange(1, ny+1, dtype=int):       
    for ix in np.arange(1 if bc_type == 0 else 0, nx, dtype=int):
        
        if scenario == 3:
            sw_in = Solar_rad(hour_angle[iy], iy)
            if ix == 0:
                #Surface temp. calc.
                ghost_cell[iy-1] = Heat_flux_surface(y[iy], ix, iy, temp, sw_in)
                temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + latent[ix, iy-1] 
            else:
                # Snowpack temp. calc.
                temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + Solar_extinction(x[ix],sw_in) + latent[ix, iy-1]  
                
        elif scenario == 4:
            temp[ix, iy] = Heat_flow(ix, iy, temp, bc_type, ghost_cell) + Solar_extinction(x[ix],SW_scaled[iy]) + latent[ix, iy-1]
              
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

net_growth = np.sum(fg, axis = 1)




############################    Data output   ############################ 

# Sensible heat calculation
if scenario == 3:
    sensible_heat = T_array[4,:]
    LW_out = T_array[5,:]
elif scenario == 4:
    sw_srf = SW_net*(1 - mt.exp(-sw_k*dx))
    sensible_heat = (((T_array[2,:]-T_array[0,:])/dx)*k) - sw_srf - LW_net
    
    LW_out = LW_net #FIX this should be changed

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
        data_titles2[0]: temp[:,0],
        }
    df = pd.DataFrame(data2)
    df.to_excel(output_file2, index=False)
    print("Wrote to file:", output_file2)
 
    
 
# Net growth to df
if run_number == 1:
    ng_data1 = {ng_title: net_growth,}
    ng_hub = pd.DataFrame(ng_data1)

else:    
    ng_data2 = {ng_title: net_growth}
    ng_df2 = pd.DataFrame(ng_data2)
    ng_hub = pd.concat([ng_hub, ng_df2], axis=1)
# Join the DataFrames by columns


temp_mean = np.mean(temp[0,2881:-1])
print('Diurnal mean surface temperature', temp_mean)

end_time = time.time()
print(f"Simulation complete. Runtime: {(end_time-start_time):.2f} seconds")


#%%
##########################     Plots    ########################## 
xticks = np.arange(0,-20,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx)




### Plot-Block: Can be run independently from model routine in Spoider  
#%%
#Surface temp
plt.figure(figsize=(9, 6))
plt.plot(y,rad_data, label= "Measured") # Plot every whole hour  
plt.plot(y, sp_temp[0,11520:], label= "Spin-up final 3-days") # Plot every whole hour   
plt.title('Surface temp')
plt.xlabel('Temperature [°C] ')
plt.ylabel('Depth [cm]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()

#%%
if scenario == 3:
    #Solar plot
    plt.plot(y,  solar_array[4,:], label= "Calc SW net") # Plot every whole hour    
    #plt.plot(y, SW_net, label= "Measured SW net") # Plot every whole hour 
    plt.plot(y, LW_out, label= "Calc LW out") # Plot every whole hour 
    plt.plot(y, sensible_heat, label= "Calc sensible heat") # Plot every whole hour
    plt.title('Shortwave rad')
    plt.xlabel('sec')
    plt.ylabel('W/m2')
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()
else:
    #Solar plot
    plt.plot(y,  solar_array[4,:], label= "Calc SW net") # Plot every whole hour    
    plt.plot(y, SW_net, label= "Measured SW net") # Plot every whole hour 
    plt.plot(y, SW_scaled, label= "Scaled SW net") # Plot every whole hour
    plt.plot(y, LW_net, label= "Calc LW out") # Plot every whole hour 
    plt.plot(y, sensible_heat, label= "Calc sensible heat") # Plot every whole hour
    plt.title('Shortwave rad')
    plt.xlabel('sec')
    plt.ylabel('W/m2')
    plt.legend()
    plt.grid(alpha=0.5)
    plt.show()


#%%
#Temperature
plt.figure(figsize=(10, 6))
plt.plot(y, temp[20,:], label= "Snow 10cm") # Plot every whole hour  
plt.plot(y, tinytag_A, label= "Tinytag A -10cm") # Plot every whole hour  
plt.plot(y,rad_data, label= "Radiometer surface temp") # Plot every whole hour  
plt.title(f'Tinytag vs snow temp. SW k={sw_k}, con 0.14 & No latent, warmer IC')
plt.xlabel('Seconds')
plt.ylabel('Temperature [°C]')
plt.legend()
plt.grid(alpha=0.5)
plt.show()


#%%
# Slide show mania

xticks = np.arange(0,-22,-2) #FIX change x-scale by MAX-MIN insted of fixed values
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

ax2.set_title("Facet growth rate")
ax2.set_xlabel("Time [h]")
ax2.set_ylabel("Facet growth rate [nm/s]")
ax2.set_yticks(np.arange(-4, 4, 0.5))
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
    line2.set_ydata(fgr[depth_idx, :])  # Update y-data (FGR at selected depth)
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
xticks = np.arange(0,-20,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx)

fig, ax = plt.subplots(3, 3, figsize=(24, 12))

# Spinup
if spin_up_has_occurred == 1:
    
    for p in np.arange(0, sp_ny+1, h*sp_pisp):
        ax[0, 0].plot(sp_temp[:, p], x, label= f"{sp_y[p]/3600} Hours")     
    ax[0, 0].set_title('Temperature profiles spin-up')
    ax[0, 0].set_xlabel('Temperature °C')
    ax[0, 0].set_ylabel('Depth [cm]')
    ax[0, 0].invert_yaxis()
    ax[0, 0].legend(fontsize=8)
    ax[0, 0].grid(alpha=0.5)
    ax[0, 0].set_xticks(xticks)
    
    
    ax[0, 1].plot(sp_temp[:,-1], x, label= f"Time {sp_y[-1]/3600} h") 
    ax[0, 1].set_title('Spin-up end state used for IC') 
    ax[0, 1].set_xlabel('Temperature [°C]')  
    ax[0, 1].set_ylabel('Depth [cm]')  
    ax[0, 1].invert_yaxis()
    ax[0, 1].legend()  
    ax[0, 1].grid(alpha=0.5)


# # Surface BC
# ax[0, 1].plot(y_hours, bc)
# ax[0, 1].set_title('Temperature surface forcing')
# ax[0, 1].set_ylabel('Temperature [°C]')
# ax[0, 1].set_xlabel('Time [h]')
# ax[0, 1].grid(alpha = 0.5)


# Temperature plot
for p in np.arange(0, ny+1, h*pisp):
    ax[0, 2].plot(temp[:pld, p], x[:pld], label=f'{y_t[p]}')
ax[0, 2].set_title('Temperature')
ax[0, 2].set_xlabel('Temperature [°C]')
ax[0, 2].set_ylabel('Depth [cm]')
ax[0, 2].invert_yaxis()
ax[0, 2].legend(fontsize=8)
ax[0, 2].grid(alpha=0.5)
ax[0, 2].set_xticks(xticks)

# Vapor Pressure plot
for p in np.arange(0, ny+1, h*pisp):
    ax[1, 0].plot(vp[:pld, p], x[:pld], label=f'{y_t[p]}')
ax[1, 0].set_title('Vapor Pressure')
ax[1, 0].set_xlabel('Vapor Pressure [mb]')
ax[1, 0].set_ylabel('Depth [cm]')
ax[1, 0].invert_yaxis()
ax[1, 0].legend()
ax[1, 0].grid(alpha=0.5)

# Vapor Pressure over time plot
ax[1, 1].plot(y_hours, vp[0, :], label='Surface vp')
ax[1, 1].set_title('Surface Vapor Pressure')
ax[1, 1].set_xlabel('Time [h]')
ax[1, 1].set_ylabel('Vapor Pressure [mb]')
ax[1, 1].legend()
ax[1, 1].grid(alpha=0.5)

# Vapor Pressure Gradient (VPG) plot
ax[1, 2].plot(y_hours, vpg[0, :], label='Surface vpg')
ax[1, 2].set_title('Vapor Pressure Gradient')
ax[1, 2].set_xlabel('Time [h]')
ax[1, 2].set_ylabel('Vapor Pressure Gradient [Pa/cm]')
ax[1, 2].legend()
ax[1, 2].grid(alpha=0.5)

# Facet Growth Rate plot
ax[2, 0].plot(y_hours, fgr[0, :], label='Surface fgr')
ax[2, 0].plot(y_hours, fgr[-1, :], label='Bottom fgr')
ax[2, 0].plot(y_hours, fgr[10, :], label='5 cm fgr')
ax[2, 0].set_title('Facet Growth Rate')
ax[2, 0].set_xlabel('Time [h]')
ax[2, 0].set_ylabel('Facet Growth Rate [nm/s]')
ax[2, 0].legend()
ax[2, 0].grid(alpha=0.5)

# Net Growth plot
ax[2, 1].plot(net_growth, x[0:-1], label='Net Growth')
ax[2, 1].set_title('Net Facet Growth')
ax[2, 1].set_xlabel('Net Growth [mm]')
ax[2, 1].set_ylabel('Depth [cm]')
ax[2, 1].invert_yaxis()
ax[2, 1].grid(alpha=0.5)

#Net Growth near surface
ax[2, 2].plot(net_growth[:pld], x[0:pld], label='Net growth near surface')
ax[2, 2].set_title('Net Facet Growth Near Surface')
ax[2, 2].set_xlabel('Net Growth [mm]')
ax[2, 2].set_ylabel('Depth [cm]')
ax[2, 2].invert_yaxis()
ax[2, 2].grid(alpha=0.5)


plt.tight_layout()
plt.show()
#%%
xticks = np.arange(0,-20,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx)
fig, ax = plt.subplots(1, 2, figsize = (12, 6), gridspec_kw={'width_ratios': [2, 1]})
#Temperature plot
cmap = plt.get_cmap("tab20")  # Alternatives: "viridis", "plasma", "tab10", "tab20", "Set3"
colors = [cmap(i/ len(np.arange(h*24, h*48, h*pisp))) for i in range(len(np.arange(h*24, h*48, h*pisp)))] 


for i, p in enumerate (np.arange(h*24, h*48, h*pisp)): # Plots day 2. FIX, make better
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax[0].plot(temp[:pld, p], x[:pld], label=time_only, color=colors[i])

ax[0].set_xlabel('Temperature [°C]')
ax[0].set_ylabel('Depth [cm]')
ax[0].invert_yaxis()
ax[0].legend(fontsize=11, title="31.03.2022")
ax[0].grid(alpha=0.5)
ax[0].set_xticks(xticks)
ax[0].xaxis.set_label_position('top')
ax[0].xaxis.tick_top()

#Net growth near surface
ax[1].plot(ng_hub['IC'][:pld], x[:pld], label='con 1.4.')
# ax[1].plot(ng_hub['IC -30%'][:pld], x[:pld], label='IC -30%')
# ax[1].plot(ng_hub['IC +30%'][:pld], x[:pld], label='IC +30%')
ax[1].legend(fontsize=11)
ax[1].set_xlabel("Net 'Facetedness' [mm]")
ax[1].set_ylabel('Depth [cm]')
ax[1].invert_yaxis()
ax[1].grid(alpha=0.5)
ax[1].xaxis.set_label_position('top')
ax[1].xaxis.tick_top()