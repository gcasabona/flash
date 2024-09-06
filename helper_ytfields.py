import yt
from yt import derived_field
import numpy as np
from yt.fields.derived_field import ValidateParameter,ValidateDataField
from yt.units import dimensions as dims
from yt.fields.local_fields import local_fields


from .helper_tabulateddata import flameprops
fn = "../../data/flameprops/flame_properties.rho.7e6-5e9.XC0.5.dat"
flameprops = flameprops(fn)

def _soundspeed(field, data):
    return np.sqrt(data["gamc"] * data["pres"] / data["dens"])
yt.add_field(("flash", "soundspeed"), function=_soundspeed, units="cm/s",sampling_type="cell")
def _density(field, data):
    return data["dens"]
yt.add_field(("flash", "density"), function=_density, units="g/cm**3",sampling_type="cell",display_name=r"$\rho$")
def _temperature(field, data):
    return data["temp"]
yt.add_field(("flash", "temperature"), function=_temperature, units="K",sampling_type="cell")
def _flamespeed(field, data):
    return data["fspd"]*yt.units.cm/yt.units.s
yt.add_field(("flash", "flamespeed"), function=_flamespeed, units="cm/s",sampling_type="cell")

def _op2mag(field, data):
    """Actually not really op2mag as multiplied by c2_h correction factor in FLASH"""
    return data["turb"]*yt.units.cm/yt.units.s
yt.add_field(("flash", "op2mag"), function=_op2mag, display_name=r"$\mathrm{mag}\left[\mathrm{OP}_2\right]$",units="cm/s",sampling_type="cell")

def _uprime(field, data):
    return data["turb"]*yt.units.cm/yt.units.s
yt.add_field(("flash", "uprime"), function=_uprime, display_name=r"$u'_\Delta$",units="cm/s",sampling_type="cell")

def _op2mag_fldtmask(field,data):
    """op_mag2 on flame surface"""
    return np.where(data[("flash","fldtmask")]>0.0,data[("flash","op2mag")],0)
yt.add_field(('flash','op2mag_fldtmask'), function=_op2mag_fldtmask,
              display_name=r"$\mathrm{mag}\left[\mathrm{OP}_2\right]$ (km/s)", sampling_type="cell",force_override=True)


def _eddytime(field,data):
    """eddy turnover time"""
    return 4*data["dx"]/data[("flash","uprime")]
yt.add_field(('flash','eddytime'), function=_eddytime,units="s",
              display_name=r"$\tau_\mathrm{eddy}=\frac{\Delta}{u'}$", sampling_type="cell",force_override=True)#

def _SToverSJ(field, data):
    fprops = data.get_field_parameter("flameprops",default=None)
    #alpha = fprops.get_value(data["density"],"density jump at the end of 12C-burning zone")
    alpha = fprops.get_value(data["density"],"12C fluid expansion factor")
    
    return data["flamespeed"]/(data["soundspeed"]/alpha)
yt.add_field(("flash", "SToverSCJ"), function=_SToverSJ, units="",sampling_type="cell",display_name=r"$S_\mathrm{T}/S_\mathrm{CJ}$",
              validators=[ValidateParameter(['flameprops'],{"flameprops":[flameprops]})])

def _velocityz(field, data):
    return data["velz"]
yt.add_field(("flash", "velocityz"), function=_velocityz, units="cm/s",sampling_type="cell",display_name="z-Velocity")
def _velocityx(field, data):
    return data["velx"]
yt.add_field(("flash", "velocityx"), function=_velocityx, units="cm/s",sampling_type="cell",display_name="x-Velocity")

def _velocityr_mag(field, data):
    return np.sqrt(data["vely"]**2+data["velz"]**2)
yt.add_field(("flash", "velocityr_mag"), function=_velocityr_mag, units="cm/s",sampling_type="cell",display_name="r-Velocity magnitude")

def _velocity_mag(field, data):
    return np.sqrt(data["velx"]**2+data["vely"]**2+data["velz"]**2)
yt.add_field(("flash", "velocity_mag"), function=_velocity_mag, units="cm/s",sampling_type="cell",display_name="r-Velocity magnitude")
def _velocity_turboverbulk(field,data):
    unitfix = yt.units.cm/yt.units.s
    return data["turb"]*unitfix/data["velocity_mag"]
yt.add_field(("flash", "velocity_turboverbulk"), function=_velocity_turboverbulk, units="",sampling_type="cell",display_name="mag[$OP_2$]/mag[vel]")

def _Umax(field,data):
    #if "fprops" not in locals():
    fprops = data.get_field_parameter("flameprops",default=None)
    rho = data[("flash","density")]
    S_L = fprops.get_value(rho,'laminar flame speed')*yt.units.cm/yt.units.s
    delta_L = fprops.get_value(rho,'width of the 12C-burning zone')*yt.units.cm
    alpha = fprops.get_value(rho,'12C fluid expansion factor')
    I_M = 1.0 # TODO: Eventually replace by better estimate.
    l = 4*data["dx"]# all cells are cubic, thus direction does not matter.
    # factor 4 due to stencil size of 2nd derivative as written in Jackson+14.
    return alpha*I_M*S_L*(l/delta_L)**(1.0/3.0)

yt.add_field(('flash','Umax'), function=_Umax, units="cm/s",
              display_name=r"$U_l^\mathrm{max}$", sampling_type="cell",force_override=True,
              validators=[ValidateParameter(['flameprops'],{"flameprops":[flameprops]})])

def _Uratio(field,data):
    """Ratio of U_l over U_l^max. For performance currently using FLASH Turb over yt-field op_mag2"""
    return data[("flash","turb")]*yt.units.cm/yt.units.s/data[("flash","Umax")]
yt.add_field(('flash','Uratio'), function=_Uratio,
              display_name=r"$U_l/U_l^\mathrm{max}$", sampling_type="cell",force_override=True)




def _fldtmask(field,data):
    """Just some simple masking for fldt"""
    thresh = 10.0
    return np.where(data[("flash","fldt")]>thresh,1.0,0.0)
yt.add_field(('flash','fldtmask'), function=_fldtmask,
              display_name=r"$\phi$ mask", sampling_type="cell",force_override=True)


def _LCJref(field,data):
    u = data.ds.units
    #u = yt.units
    fprops = data.get_field_parameter("flameprops",default=None)
    rho = data[("flash","density")]
    S_L = fprops.get_value(rho,'laminar flame speed')*u.cm/u.s
    cs = data[("flash","soundspeed")]
    alpha = fprops.get_value(rho,'12C fluid expansion factor')
    I_M = 1.0 # TODO: Eventually replace by better estimate.
    l = 4.0*data["dx"]# all cells are cubic, thus direction does not matter.
    # factor 4 due to stencil size of 2nd derivative as written in Jackson+14.
    U_l = data[("flash","turb")]*u.cm/u.s
    #result = (alpha*I_M*S_L)**2*cs*(l/U_l**3)
    a = (alpha*I_M*S_L)**2
    b = cs*(l/U_l**3)
    result = a*b
    return result 

yt.add_field(('flash','LCJref'), function=_LCJref, #units="cm",
              display_name=r"$L_\mathrm{CJ}^\mathrm{ref}$", sampling_type="cell",force_override=True,
              dimensions=dims.length,units="auto",
              validators=[ValidateParameter(['flameprops'],{"flameprops":[flameprops]})])

def _LCJmin(field,data):
    u = data.ds.units
    fprops = data.get_field_parameter("flameprops",default=None)
    rho = data[("flash","density")]
    S_L = fprops.get_value(rho,'laminar flame speed')*u.cm/u.s
    cs = data[("flash","soundspeed")]
    alpha = fprops.get_value(rho,'12C fluid expansion factor')
    I_M = 1.0 # TODO: Eventually replace by better estimate.
    l = 4.0*data["dx"]# all cells are cubic, thus direction does not matter.
    # factor 4 due to stencil size of 2nd derivative as written in Jackson+14.
    U_l = data[("flash","turb")]*u.cm/u.s
    delta_L = fprops.get_value(rho,'width of the 12C-burning zone')*u.cm
    result = (delta_L*cs)/(alpha*I_M*S_L)
    return result

yt.add_field(('flash','LCJmin'), function=_LCJmin, #units="cm",
              display_name=r"$L_\mathrm{CJ}^\mathrm{min}$", sampling_type="cell",force_override=True,
              dimensions=dims.length,units="auto",
              validators=[ValidateParameter(['flameprops'],{"flameprops":[flameprops]})])


def _LCJ(field,data):
    return np.maximum(data[("flash","LCJref")],data[("flash","LCJmin")])
yt.add_field(('flash','LCJ'), function=_LCJ, dimensions=dims.length,units="auto",
              display_name=r"$L_\mathrm{CJ}$", sampling_type="cell",force_override=True)

def _Nddt(field,data):
    return (4*data["dx"]/data[("flash","LCJ")])**(11.0/3.0)
yt.add_field(('flash','Nddt'), function=_Nddt,
              display_name=r"$N_\mathrm{DDT}$", sampling_type="cell",force_override=True)


def _Lratio_LCJ(field,data):
    """Ratio of L over LCJ. For performance currently using FLASH Turb over yt-field op_mag2"""
    return 4*data["dx"]/data[("flash","LCJ")]
yt.add_field(('flash','Lratio_LCJ'), function=_Lratio_LCJ,
              display_name=r"$L/L_\mathrm{CJ}$", sampling_type="cell",force_override=True)
def _Lratio_LCJ_fldtmask(field,data):
    """Ratio of L over LCJ. For performance currently using FLASH Turb over yt-field op_mag2"""
    return np.where(data[("flash","fldtmask")]>0.0,data[("flash","Lratio_LCJ")],0)
yt.add_field(('flash','Lratio_LCJ_fldtmask'), function=_Lratio_LCJ_fldtmask,
              display_name=r"$L/L_\mathrm{CJ}$", sampling_type="cell",force_override=True)

def _Lratio(field,data):
    """Ratio of L over max(LCJref,LCJmin). For performance currently using FLASH Turb over yt-field op_mag2"""
    return 4*data["dx"]/data[("flash","LCJ")]
yt.add_field(('flash','Lratio'), function=_Lratio_LCJ,
              display_name=r"$L/\max\left(L_\mathrm{CJ}^\mathrm{min},L_\mathrm{CJ}^\mathrm{ref}\right)$", sampling_type="cell",force_override=True)
#def _Lratio_fldtmask(field,data):
#    """Ratio of L over LCJ. For performance currently using FLASH Turb over yt-field op_mag2"""
#    return np.where(data[("flash","fldtmask")]>0.0,data[("flash","Lratio")],0)
#yt.add_field(('flash','Lratio_fldtmask'), function=_Lratio_fldtmask,
#              display_name=r"$L/\max\left(L_\mathrm{CJ}^\mathrm{min},L_\mathrm{CJ}^\mathrm{ref}\right)$", sampling_type="cell",force_override=True)

def _LratioCJmin(field,data):
    """Ratio of LCJ over LCJmin. For performance currently using FLASH Turb over yt-field op_mag2"""
    return data[("flash","LCJ")]/data[("flash","LCJmin")]
yt.add_field(('flash','LratioCJmin'), function=_LratioCJmin,
              display_name=r"$L_\mathrm{CJ}/L_\mathrm{CJ}^\mathrm{min}$", sampling_type="cell",force_override=True)




## !!!!!!!!!!!!!!!!!!
## lots of junk below
## !!!!!!!!!!!!!!!!!!
def _velzoversoundspeed(field, data):
    return data["velocityz"]/data["soundspeed"]
yt.add_field(("flash", "velzoversoundspeed"), function=_velzoversoundspeed, units="",
              sampling_type="cell",display_name=r"$v_\mathrm{z}/c_s$")



def _velzoverSJ(field, data):
    return data["velocityz"]/(data["soundspeed"]/1.1)
yt.add_field(("flash", "velzoverSJ"), function=_velzoverSJ, units="",sampling_type="cell")


def _velzoverSJabs(field, data):
    return np.abs(data["velocityz"]/(data["soundspeed"]/1.1))
yt.add_field(("flash", "velzoverSJabs"), function=_velzoverSJabs, units="",sampling_type="cell")

def _TURBoverSJ(field, data):
    return data["turb"]/(data["soundspeed"]/1.1)
yt.add_field(("flash", "TURBoverSCJ"), function=_TURBoverSJ, units="",sampling_type="cell")
## !!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!
## !!!!!!!!!!!!!!!!!!

### certain fields, we only want to evaluate on the flame surface.
flamesurfacefields = ["op2mag","Nddt","Lratio"]
for f in flamesurfacefields:
    def fnc(field,data):
        mask = data[("flash","fldtmask")]>0.0
        return np.where(mask,data["flash",f],0.0)
    units = local_fields[("flash",f)].units
    display_name = local_fields[("flash",f)].display_name
    yt.add_field(('flash',f+'_fldtmask'), function=fnc,units=units,display_name=display_name,
                   sampling_type="cell",force_override=True)


# Summarize Plot Parameters
fieldplotprops = {}
fieldplotprops[("flash","soundspeed")] = dict(cmap="cividis",cmap_set_over="white",log=True,zlim=[5e8,1e9])
fieldplotprops[("flash","density")] = dict(cmap="Blues",cmap_set_under="white",log=True,zlim=[1e3,2e9])
fieldplotprops[("flash","temperature")] = dict(cmap="magma",log=True,zlim=[1e6,1e10])
fieldplotprops[("flash","velocityz")] = dict(cmap="coolwarm",log=True,zlim=[-2e9,2e9],zlinthresh=1e7)
fieldplotprops[("flash","velocityx")] = dict(cmap="coolwarm",log=True,zlim=[-2e8,2e8],zlinthresh=1e7)
fieldplotprops[("flash","velocityr_mag")] = dict(cmap="magma",log=True,zlim=[1e4,1e9])
fieldplotprops[("flash","velocity_mag")] = dict(cmap="magma",log=True,zlim=[1e4,1e9])
#fieldplotprops[("flash","SToverSCJ")] = dict(cmap="bwr",log=False,zlim=[0,2])
fieldplotprops[("flash","SToverSCJ")] = dict(cmap="Reds",log=True,zlim=[1e-3,1.0],cmap_set_over="blue")
fieldplotprops[("flash","velzoversoundspeed")] = dict(cmap="seismic",log=False,zlim=[-1,1],cmap_set_over="green",cmap_set_under="green") #,zlinthresh=0.01,zlim=[-2,2]
fieldplotprops[("flash","turb")] = dict(cmap="viridis",log=True,zlim=[1e6,2e8])
fieldplotprops[("flash","op2mag")] = dict(cmap="viridis",log=True,zlim=[1e6,2e8])
fieldplotprops[("flash","Lratio")] = dict(cmap="plasma",log=True,zlim=[1e-2,1e4])


for f in flamesurfacefields:
    print("add "+f)
    if ("flash",f) in fieldplotprops:
        fieldplotprops[("flash",f+"_fldtmask")] = fieldplotprops[("flash",f)]

# make flash category of field also default for plot parameter dictionary
for t in list(fieldplotprops.keys()):
    if type(t)==tuple:
        if t[0]=="flash":
            fieldplotprops[t[1]] = fieldplotprops[t]

