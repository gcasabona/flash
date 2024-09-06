import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import pandas as pd
import os
import yt

%run -i python/helper_tabulateddata.py
fn = "data/flameprops/flame_properties.rho.7e6-5e9.XC0.5.dat"
fprops = flameprops(fn)

#Function that computes for sound speed
def _CS(field, data):
    pressure = data[('flash','pres')]
    density = data[('flash','dens')]
    gamma = data[('flash','gamc')]
    
    return (gamma*pressure/density)**(1./2)

#Function that computes for the Chapman-Jouguet length
def _LCJ(field, data):
    dens = data[('flash','dens')] #density
    spd_sound = data[('flash','sndspd')] #speed of sound
    S_L = fprops.get_value(dens,'laminar flame speed')*yt.units.cm/yt.units.s #laminar flame speed
    alpha = fprops.get_value(dens,'12C fluid expansion factor') #fuel-to-ash burn ratio
    l = data["dx"]
    U_l = data[('flash','turb')]*yt.units.cm/yt.units.s #turbulence intensity
    I_M = 1.0 #coefficient of turbulent stretch, for thermonuclear burning =1
    
    return (alpha*I_M*S_L)**2*spd_sound*(l/U_l**3)

#Function that computes for the ratio between the integral scale and the Chapman-Jouguet Length

def _Lratio1(field,data):
    return data["dx"]/data[('flash','LCJ')]

#Function that computes for the ratio between the Chapman-Jouguet Length and the laminar flame thickness

def _Lratio2(field,data):
    dens = data[('flash','dens')]
    delta_L = fprops.get_value(dens,'width of the 12C-burning zone')*yt.units.cm # laminar flame thickness
    L_CJ = data[('flash','LCJ')]
    
    return (L_CJ/delta_L)

#Function that computes for the ratio between turbulent velocity and the turbulent velocity 
#that provides maximum burning rate which could lead to tDDT

def _Uratio(field,data):
    dens = data[('flash','dens')] #density
    spd_sound = data[('flash','sndspd')] #speed of sound
    S_L = fprops.get_value(dens,'laminar flame speed')*yt.units.cm/yt.units.s #laminar flame speed
    delta_L = fprops.get_value(dens,'width of the 12C-burning zone')*yt.units.cm # laminar flame thickness
    alpha = fprops.get_value(dens,'12C fluid expansion factor') #fuel-to-ash burn ratio
    l = 4*data["dx"]
    U_l = data[('flash','turb')]*yt.units.cm/yt.units.s #turbulence intensity
    I_M = 1.0 #coefficient of turbulent stretch, for thermonuclear burning =1
    
    return U_l/(alpha*I_M*S_L*(l/delta_L)**(1.0/3.0))

#Function that computes for the ratio between the flame speed and the Chapman-Jouguet flame speed

def _Sratio(field,data):
    dens = data[('flash','dens')] #density
    alpha = fprops.get_value(dens,'12C fluid expansion factor') #fuel-to-ash burn ratio
    
    return (data[('flash','fspd')]*yt.units.cm/yt.units.s)/(data['flash','sndspd']/alpha)

#Function that computes for the RMS velocity u_rms sqrt(1/3(u_x**2 + u_y**2 + u_z**2))

def _RMS(field,data):
    vel_x = data[('all','particle_velx')] # x-velocity
    vel_y = data[('all','particle_vely')] # y-velocity
    vel_z = data[('all','particle_velz')] # z-velocity

    return np.sqrt((1./3)*(vel_x**2 + vel_y**2 + vel_z**2))

# Add derived yt fields
yt.add_field(('flash','Sratio'), function=_Sratio,
              display_name=r"$S_T/S_\mathrm{CJ}$", sampling_type="cell",force_override=True)
yt.add_field(('flash','Uratio'), function=_Uratio,
              display_name=r"$U_l/U_l^\mathrm{max}$", sampling_type="cell",force_override=True)
yt.add_field(('flash','Lratio1'), function=_Lratio1,
              display_name=r"$L/L_\mathrm{CJ}$", sampling_type="cell",force_override=True)
yt.add_field(('flash','Lratio2'), function=_Lratio2,
              display_name=r"$L_\mathrm{CJ}/\delta_L$", sampling_type="cell",force_override=True)
yt.add_field(('flash','sndspd'),function = _CS, units="cm/s", 
             display_name = "sound speed", force_override=True)
yt.add_field(('flash','LCJ'), function=_LCJ, units="cm", 
             display_name=r"$L_\mathrm{CJ}$", sampling_type="cell",force_override=True)

#Function that produces a plot panel of global observables through time
def MakeGlobal(file):
    data_head = np.array(pd.read_csv(file,  delim_whitespace=True).columns)
    data = np.genfromtxt(file)
    data_new = np.ndarray.transpose(np.zeros(shape = np.shape(data)))
    
    for i in range(np.shape(data_new)[0]):
        for j in range(np.shape(data_new)[1]):
            data_new[i,j] = data[j,i]
            
    no = int(np.shape(data_new)[0]/3)+1
    fig, axes = plt.subplots(int(np.shape(data_new)[0]/3), 3, figsize=(15,3*no), tight_layout=True)
    axes = axes.flatten()
    
    for j in range(np.shape(data_new)[0]-1):
        axes[j].plot(data_new[0],data_new[j+1])
        axes[j].set_title(data_head[j+1])
    
    plt.savefig(os.path.splitext(file)[0]+".jpg")
    
# Function that makes multipanel plot
def MakePanel(data):
    fig = plt.figure()
    grid = AxesGrid(fig, (0.1,0.1,0.85,0.85),
                nrows_ncols = (2, 2),
                axes_pad = 1.0,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="5%",
                cbar_pad="0%")
    fields = [('flash','Sratio'),('flash','Lratio1'),('flash','Lratio2'),('flash','Uratio')]
    plot = yt.SlicePlot(data, 'x', fields,width = (3000,"km"))
    
    for i, field in enumerate(fields):
        p = plot.plots[field]
        p.figure = fig
        p.axes = grid[i].axes
        p.cax = grid.cbar_axes[i]
    
    plot._setup_plots()
    plt.savefig('multi_panel.png')

# Function that displays basic volume render of yt fields
def MakeVolume(data,field):
    im, sc = yt.volume_render(data, field, fname = 'render1.png')
    sc.camera.width = (500, 'km')
    sc.camera.switch_orientation()
    sc.save('render2.png')
    sc.show()

