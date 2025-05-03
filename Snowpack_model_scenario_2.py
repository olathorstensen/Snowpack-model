### 1D snowpack temperature simulator ###
# v1.4 
#@author: Ola Thorstensen and Thor Parmentier
# Version update:
#   - Figure for paper plot
 
    
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.widgets import Slider
from datetime import datetime, timedelta
import time
start_time = time.time()

############################    Parameters   ############################ 
runtime = 24             # Hours
dt = 30                  # Time step [seconds] (Must be a divisor of 3600)
depth = 1                # Snow depth from surface [m]
dx = 0.005               # Dist. interval [m]
pisp = 2                 # USE INTEGERS. Give hourly plot rate
b_bc = 0                # Bottom boundary condition, fixed [°C]
spin_up = 1              # [0] No spin-up, [1] Run spin-up
sp_runtime = 24*3        # Spin-up run time [Hours]
plot_depth = 0.40        # Depth shown in plots [m]
ng_to_file = 1           # Writes net_growth to file


############################    Constants   ############################   
k = 0.1439               # Thermal conductivity snow [W/m K]
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
    
    
#########################    Model domain    #########################

    # Model grid  
nx = int(depth/dx)
ny = int((runtime * 60 * 60) / dt)

temp = np.zeros([nx+1, ny+1], dtype=float)  # Main snowpack grid
vp = np.zeros([nx+1, ny+1], dtype=float)    # Vapor pressure grid
vpg = np.zeros([nx, ny+1], dtype=float)    # Vapor pressure gradient grid
fgr = np.zeros([nx, ny+1], dtype=float)    # Facet growth rate grid
fg =  np.zeros([nx, ny], dtype=float)    # Facet growth grid

    # Axis and time
x = np.linspace(0, depth*100, nx+1) # Depth axis
x_stag = x[:-1]+dx*100/2            # Staggered x axis for vpg, fgr, fg, and ng
y = np.round(np.arange(0, ny+1, 1) *dt) #/3600) # Time axis in sec (divide by 3600 for h)
y_sec = y.astype(float)
y_hours = np.round(np.arange(0, ny+1, 1) *dt)/3600

base_time = datetime.strptime("00:00 30/03/2022", "%H:%M %d/%m/%Y") 
y_t = [(base_time + timedelta(seconds=seconds)).strftime("%H:%M %d/%m/%Y") for seconds in y_sec]


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



temp[0,:] = bc
temp[-1,:] = b_bc * np.ones(ny+1, dtype=float)  #Fixed bottom bc

# Could move forcing param to "Parameter section" 
    
        
# Prints of shapes and numbers (all are numbers and none are shapes)
print('r_number (should be less than 0.5):', r)
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
        
for iy in np.arange(0, ny, dtype=int):
    for ix in np.arange(0, nx, dtype=int):
        fg[ix, iy] = Facet_growth(ix, iy)

if b_bc<0:
    net_growth_cold_bc = np.sum(fg, axis = 1)
else:
    net_growth = np.sum(fg, axis = 1)

############################    Data output   ############################ 

current_dir = os.path.dirname(os.path.abspath(__file__)) 
if ng_to_file == 1:
    output_file = os.path.join(current_dir, 'Net_growth_data.xlsx')
    data_titles = ["Scenario 2"] #f"Net growth after {sp_y[-1]} hours"
    data = {
        data_titles[0]: net_growth,
        }
    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False)
    print("Wrote to file:", output_file)

#########################################################################

end_time = time.time()
print(f"Simulation complete. Runtime: {(end_time-start_time):.2f} seconds")

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

### Big fig
fig, ax = plt.subplots(3, 3, figsize=(24, 12))

# Spinup
if spin_up == 1:
    ax[0, 0].plot(sp_temp[:, -1], x, label=f"Time {sp_y[-1]} h")
    ax[0, 0].set_title('Spin-up end state used for IC') 
    ax[0, 0].set_xlabel('Temperature [°C]')  
    ax[0, 0].set_ylabel('Depth [cm]')  
    ax[0, 0].invert_yaxis()
    ax[0, 0].legend()  
    ax[0, 0].grid(alpha=0.5)

# Surface BC
ax[0, 1].plot(y_hours, bc)
ax[0, 1].set_title('Temperature surface forcing')
ax[0, 1].set_ylabel('Temperature [°C]')
ax[0, 1].set_xlabel('Time [h]')
ax[0, 1].grid(alpha = 0.5)

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
### Figure for the paper
plt.rcParams.update({'font.size': 24})
xticks = np.arange(0,-22,-2) #FIX change x-scale by MAX-MIN insted of fixed values
pld = int(plot_depth/dx) +1
time_steps = ny + 1

#fig3, ax3 = plt.subplots(1, 2, figsize = (12, 6), gridspec_kw={'width_ratios': [2, 1]})
fig3, ax3 = plt.subplots(1, 2, figsize=(12, 6), gridspec_kw={'width_ratios': [2, 1], 'wspace': 0.1})
# Temperature plot
cmap = plt.get_cmap("tab20")  # Alternatives: "viridis", "plasma", "tab10", "tab20", "Set3"
colors = [cmap(i / len(np.arange(0, ny, h*pisp))) for i in range(len(np.arange(0, ny+1, h*pisp)))]


for i, p in enumerate(np.arange(0, ny, h*pisp)):
    time_only = (base_time + timedelta(seconds=y_sec[p])).strftime("%H:%M")
    ax3[0].plot(temp[:pld, p], x[:pld], label=f'{time_only}', color=colors[i], lw=2.5)
ax3[0].set_xlabel('Temperature [°C]')#, fontsize = 14)
ax3[0].set_ylabel('Depth [cm]')#, fontsize = 14)
ax3[0].yaxis.tick_right()  # Move y-tick labels to the right
ax3[0].yaxis.set_label_position("right")  # Move y-axis label to the right
ax3[0].set_ylabel('Depth [cm]', rotation=270, labelpad=25)
ax3[0].invert_yaxis()
ax3[0].legend()#fontsize=13)
ax3[0].grid(alpha=0.5)
ax3[0].set_xticks(xticks)
ax3[0].xaxis.set_label_position('top')
ax3[0].xaxis.tick_top()  
ax3[0].text(0.96, 0.05, "a", 
           #fontsize=15, 
           #fontweight='bold', 
           transform=ax3[0].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))

pld = int(plot_depth/dx)
# Net growth near surface
ax3[1].plot(net_growth[:pld], x_stag[0:(pld)],linestyle='dotted', color='C3', lw=5, label='0°C')
ax3[1].plot(net_growth_cold_bc[:pld], x_stag[0:(pld)],linestyle= '--', color='C0', lw=2.7, label='-10°C')
ax3[1].set_xlabel("Net 'facetedness' [mm]")#, fontsize = 14)
#ax3[1].set_ylabel('Depth [cm]')#, fontsize = 14)
#ax3[1].legend(loc='lower left')
ax3[1].set_yticklabels([])
ax3[1].invert_yaxis()
ax3[1].grid(alpha=0.5)
ax3[1].xaxis.set_label_position('top')
ax3[1].xaxis.tick_top()
ax3[1].set_xticks([-0.02, -0.01, 0, 0.01, 0.02, 0.03])
ax3[1].text(0.91, 0.05, "b", 
           #fontsize=15, 
           #fontweight='bold', 
           transform=ax3[1].transAxes,
           verticalalignment='top', 
           horizontalalignment='left', 
           bbox=dict(facecolor='white', edgecolor='black', boxstyle='square,pad=0.3'))
plt.tight_layout()

