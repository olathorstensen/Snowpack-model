### 1D snowpack temperature simulator ###
# v1.8 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
#   - Plotting post review
#   - Code cleanup
 
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.widgets import Slider
from datetime import datetime, timedelta
import time
start_time = time.time()

############################    Parameters   ############################ 
runs = 2                 # Number of model runs (spinup + main model) must be 
runtime = 24             # Hours
dt = 30                  # Time step [seconds] (Must be a divisor of 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.005               # Dist. interval [m]
pisp = 2                 # USE INTEGERS. Give hourly plot rate
b_bc = 0                 # Bottom boundary condition, fixed [°C]
spin_up = 1              # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*14        # Spin-up run time [Hours]
plot_depth = 0.40        # Depth shown in plots [m]
ng_to_file = 1           # Writes net_growth to file


############################    Snow Parameters 2nd run   ############################ 
b_bc_2 = -10                 # Bottom boundary condition, fixed [°C]


############################    Constants   ############################   
k = 0.1439               # Thermal conductivity snow [W/m K]  0.1439 
rho = 245                # density [kg/m3]
cp = 2090                # Specific heat capacity of ice [J/kg °C]
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
    fg = (fgr[ix, iy] + fgr[ix, iy + 1]) / 2 * dt / 10**6
    
    return fg



    
    

#########################    Model loop    #########################
for run in range (1, runs + 1):


    #########################    Model domain    #########################        
        # Model grid  
    nx = int(depth/dx)
    ny = int((runtime * 60 * 60) / dt)
    
    temp = np.zeros([nx+1, ny+1], dtype=float)  # Main snowpack grid
    vp = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure grid
    vpg = np.zeros([nx, ny+1], dtype=float)    # Vapor pressure gradient grid
    fgr = np.zeros([nx, ny+1], dtype=float)    # Facet growth rate grid
    fg =  np.zeros([nx, ny], dtype=float)    # Facet growth grid
    dtvpge = np.zeros_like(vpg, dtype=float)
    
        # Axis and time
    x = np.linspace(0, depth * 100, nx + 1) # Depth axis
    x_stag = x[:-1]+dx*100/2            # Staggered x axis for vpg, fgr, fg, and ng
    y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
    y_sec = y.astype(float)
    y_hours = np.round(np.arange(0, ny+1, 1) *dt)/3600
    
    base_time = datetime.strptime("00:00 30/03/2022", "%H:%M %d/%m/%Y") 
    y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M %d/%m/%Y") for seconds in y_sec]
  
    
        # Initial condition
    # Linear ic
    if run == 2:
        b_bc = b_bc_2
    ic = np.linspace(-10, b_bc, nx+1)
    temp[:, 0] = np.round(ic, 4)
    
    # Fixed bc, diurnal oscillation:
    xx = np.linspace(-90, 270, int(h * 24) + 1)
    bc = 9*np.sin(np.deg2rad(xx))-10      # Sinusoidal surface bc
    bc_dummy = np.delete(bc, 0)
    sp_bc = bc
    
    # Extends or crops 24h temperature swings to runtime length
    if runtime > 24:
        for i in np.arange(1, int(runtime/24)):
            bc = np.concatenate((bc, bc_dummy,), axis=0)
        # Add remainder
        if runtime%24 != 0:
            bc_ext = bc_dummy[0:(int(runtime%24*h))]
            bc = np.concatenate((bc, bc_ext))
    
    if runtime < 24:
        bc = bc[0:int(runtime*h)+1]  # Crops temp forcing
    
    temp[0,:] = bc
    temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)  # Fixed bottom bc
    
    print('r_number (should be less than 0.5):', r)


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
        sp_y = np.round(np.arange(0, sp_ny+1, 1) * dt/3600 , 2) #Spin-up y-axis
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
    
    
    #####################     Main Model Calculations    ##################### 
    
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
            
    for iy in np.arange(0, ny, dtype=int):
        for ix in np.arange(0, nx, dtype=int):
            fg[ix, iy] = Facet_growth(ix, iy)
            
            
            
    # DTVPG calculations
    dtvpge = np.where(vpg < -5, vpg + 5, np.where(vpg > 5, vpg - 5, 0))
    
    if run == 2:
        dtvpge_mean_cold_bc = np.mean(dtvpge, axis = 1)
    else:
        dtvpge_mean = np.mean(dtvpge, axis = 1)
        
 

    
    
    ############################    Data output   ############################ 
    current_dir = os.path.dirname(os.path.abspath(__file__)) 
    if ng_to_file == 1 and run == 1:
        output_file = os.path.join(current_dir, 'DTVPGE_data.xlsx')
        data_titles = ["Scenario 2"] #f"DTVPGE after {sp_y[-1]} hours"
        data = {
            data_titles[0]: dtvpge_mean,
            }
        df = pd.DataFrame(data)
        df.to_excel(output_file, index=False)
        print("Wrote to file:", output_file)
        

#########################################################################

end_time = time.time()
print(f"Simulation complete. Runtime: {(end_time-start_time):.2f} seconds")


#%%
### Figure for the paper
plt.rcParams.update({'font.size': 24})
pld = int(plot_depth/dx) + 1
time_steps = ny + 1


fig3, ax3 = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.1})

# Temperature plot
cmap = plt.get_cmap("tab20")  # Alternatives: "viridis", "plasma", "tab10", "tab20", "Set3"
colors = [cmap(i / len(np.arange(0, ny, h*pisp))) for i in range(len(np.arange(0, ny+1, h*pisp)))]


for i, p in enumerate(np.arange(0, ny, h*pisp)):
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax3[0].plot(temp[:pld, p], x[:pld], label=f'{time_only}', color=colors[i], lw=2.5)
ax3[0].set_xlabel('Temperature [°C]')
ax3[0].set_ylabel('Depth [cm]')
ax3[0].yaxis.tick_right() 
ax3[0].yaxis.set_label_position("right")  
ax3[0].set_ylabel('Depth [cm]', rotation=270, labelpad=25)
ax3[0].invert_yaxis()
ax3[0].legend()
ax3[0].grid(alpha=0.5)
ax3[0].set_xticks(np.arange(0,-22,-2))
ax3[0].xaxis.set_label_position('top')
ax3[0].xaxis.tick_top()  
ax3[0].text(0.968, 0.05, "a",  
           transform=ax3[0].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

pld = int(plot_depth/dx) + 1

# DTVPGE near surface
ax3[1].plot(dtvpge_mean[:pld], x_stag[0:(pld)],linestyle='dotted', color='C3', lw=5, label='0°C')
ax3[1].plot(dtvpge_mean_cold_bc[:pld], x_stag[0:(pld)],linestyle= '--', color='C0', lw=2.7, label='-10°C')
ax3[1].set_xlabel("DTVPGE [Pa/cm]")
ax3[1].legend(loc='lower right', bbox_to_anchor=(0.9, 0))
ax3[1].set_yticklabels([])
ax3[1].invert_yaxis()
ax3[1].grid(alpha=0.5)
ax3[1].xaxis.set_label_position('top')
ax3[1].xaxis.tick_top()
ax3[1].text(0.935, 0.05, "b", 
            transform=ax3[1].transAxes,
            verticalalignment='top', 
            horizontalalignment='left', 
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
plt.tight_layout()

#%%
# Plot of results:
#   %matplotlib qt                # To be pasted in console before plotting
    
    
xticks = np.arange(0,-20,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx) + 1
time_steps = ny + 1

### Sliders
# Temp
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

# Create slider axis
slider_ax1 = plt.axes([0.2, 0.05, 0.6, 0.03])  # [left, bottom, width, height]
time_slider = Slider(slider_ax1, "Time Step", 0, time_steps - 1, valinit=0, valstep=1)

def update_time(val):
    time_idx = int(time_slider.val)
    line1.set_xdata(temp[:pld, time_idx])
    ax1.legend([f'Time: {y_t[time_idx]}'], loc='upper right')
    # fig1.canvas.draw_idle()  # Redraw the figure. Might not actually need this.

time_slider.on_changed(update_time)
plt.show()


# Growth rate
fig2, ax2 = plt.subplots(figsize=(9, 6))
plt.subplots_adjust(bottom=0.2)  # Space for slider
line2, = ax2.plot(y_hours, fgr[0, :], label = f'Depth: {x[0]}')

ax2.set_title("Facet growth rate")
ax2.set_xlabel("Time [h]")
ax2.set_ylabel("Facet growth rate [nm/s]")
ax2.legend()
ax2.grid(alpha=0.5)

slider_ax2 = plt.axes([0.2, 0.05, 0.6, 0.03])  # [left, bottom, width, height]
depth_slider = Slider(slider_ax2, "Depth", 0, len(x) - 2, valinit=0, valstep=1)

def update_depth(val):
    depth_idx = int(depth_slider.val)
    line2.set_ydata(fgr[depth_idx, :])  # Update y-data (FGR at selected depth)
    ax2.legend([f'Depth: {x[depth_idx]:.1f} cm'], loc='upper right')  
    fig2.canvas.draw_idle()  # Redraw the figure. Might not need this

depth_slider.on_changed(update_depth)
plt.show()

