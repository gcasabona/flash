# estimate runaway temperature
import math
import numpy as np
import matplotlib.pyplot as plt

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

rho = 1.e5
rho9 = rho / 1.e9
#rho9 = 1.e-3 # 1.e-2
rholist = [1.e4, 1.e5, 1.e6]
rhoarr = np.asarray (rholist)

x4 = 1.   #x12 = 0.25 # 0.5
v0 = 1.e8 #1.78e5  # 1.e3 km/s
L = 1.e7   # 100 km



#vr = v0 * (r / L)**(1./3.)
#tedd = r / vr

Tcrit = 1.6e9 # Temperature where enuc = eturb
Tmin  = 1.e8 #4.05e8  # Minimum temperature on plot
Tmax =  4.e9 # 3.6e9 # Maximum temperature on plot
#r = 1.e5   # 1 km
#v = v0 * (r / L)**(1./3.)

# Peak temperatures of C/O mergers
# M. Dan, S. Rosswog, M. Bruggen, P. Podsiadlowski, 2014
# MNRAS, 438, 14 (arXiv:1308.1667)
# 3.229e8 - 1.0407e9
vertline = np.logspace (-20, 0, num = 100)
temp = np.logspace (math.log10 (Tmin), math.log10 (Tmax), num = 100)
t9   = temp / 1.e9
t8   = temp / 1.e8

horizontal = temp / temp

ax = plt.axes()
ax.arrow(0, 0, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')

O = 9.37e+09/t9**(2./3.)*np.exp(-39.757/t9**(1./3.)-(t9/1.586)**2)+6.21e+01/t9**(3./2.)*np.exp(-10.297/t9)+5.38e+02/t9**(3./2.)*np.exp(-12.226/t9)+1.30e+01*t9**2*np.exp(-20.093/t9)

C = 1.04e+08/(t9**2.)/((1.00+0.0489/t9**(2./3.))**2.)*np.exp(-32.120/t9**(1./3.) - (t9/3.496)**2.) + 1.76e+08/(t9**2.)/((1.00+0.2654/t9**(2./3.))**2.)*np.exp(-32.120/t9**(1./3.)) + 1.25e+03/(t9**(3./2.))*np.exp(-27.499/t9) + 1.43e-02*(t9**5.)*np.exp(-15.541/t9)  

Be = 1.30e+02/t9**(3./2.)*np.exp(-3.3364/t9)+2.510e+07/t9**(2./3.)*np.exp(-23.570/t9**(1./3.)-(t9/0.235)**2)*(1.0+0.018*t9**(1./3.)+5.249*t9**(2./3.)+0.650*t9+19.176*t9**(4./3.)+6.034*t9**(5./3.)) 

He = 7.40e+05/t9**(3./2.)*np.exp(-1.0663/t9)+4.164e+09/t9**(2./3.)*np.exp(-13.490/t9**(1./3.)-(t9/0.098)**2)*(1.0+0.031*t9**(1./3.)+8.009*t9**(2./3.)+1.732*t9+49.883*t9**(4./3.)+27.426*t9**(5./3.))

#enuc =  3.9957616278911012E+018 * rho9**(2.5) * (x12**2.) * (O + Be + C + He)
#enuc_5 = 3.9957616278911012E+018 * (rho9*1.e-1)**(2.5) * (x12**2.) * (O + Be + C + He)

#enuc = 3.e11 * x4**3. * rho9**2. * t8**(-3.) * np.exp (-43.5 /t8) #  5.1e8 * rho9**(2.) * (x4**2.) * (t9**(-3.)) * np.exp(-4.4027/t9) # (O + Be + C + He) * (t9**(-3.)) * np.exp(-4.4027/t9)
#enuc_5 = 3.e11 * x4**3. * (rho9*1.e-1)**2. * t8**(-3.) * np.exp (-43.5 /t8) # 5.1e8 * (rho9*1.e-1)**(2.) * (x4**2.) * (t9**(-3.)) * np.exp(-4.4027/t9) # (O + Be + C + He) * (t9**(-3.)) * np.exp(-4.4027/t9)
#2.06e46

eturb = v0**3. / L

#ratio = enuc / eturb
#ratio_5 = enuc_5 / eturb

for rho in rhoarr :
   enuc = 3.9e11 * rho**2 * x4**3 * t8**(-3) * np.exp(2.76e-3 * rho**(1./2.) * t8**(-3./2.)) * np.exp(-42.94/t8)  #3.e11 * x4**3. * rho**2. * t8**(-3.) * np.exp (-43.5 /t8)
   ratio = enuc / eturb
   plt.loglog (temp, ratio, 'k', linewidth = 2.5)

#plt.loglog (temp, ratio, 'k', linewidth = 2.5)
#plt.loglog (temp, ratio_5, 'k', linewidth = 2.5)
#plt.fill_between (temp, vertline,  alpha = 0.5, hatch = '/')
plt.loglog (temp, horizontal, '--', linewidth = 1.5)
plt.xlabel ("Temperature (K)")
#plt.ylabel ("Ratio $\epsilon_{\rm nuc}/ \epsilon_{\rm turb}$")
plt.ylabel ("Ratio Specific Nuclear to Turbulent Dissipation Rates")
#plt.axvline(x = 2.2e9)
#plt.axvline(x=1.0407e9)
#plt.axvline(x=3.229e8)
#plt.axvspan(3.229e8, 1.0407e9, alpha=0.25, color='dimgrey')
plt.text (Tmin, 3.e0, "$\epsilon_{turb} = \epsilon_{nuc}$")
plt.text (Tcrit, 1.e-5, "Nuclear-Dominated")
ax.annotate("", xy=(Tcrit, 1.e-6), xytext=(3.e9, 1.e-6),
    arrowprops=dict(arrowstyle="<-"))
ax.annotate("", xy=(Tcrit, 2.e-13), xytext=(1.e9, 2.e-13),
    arrowprops=dict(arrowstyle="<-"))
plt.text (Tcrit, 1.e-12, "Turbulent-Dominated", ha = "right")
plt.text (3.0e8, 1.e-15, "WD Merger Peak Temp")
plt.axvspan(Tcrit, Tmax, alpha=0.25, color='red')
ax = plt.axes()
ax.annotate("", xy=(3.22e8, 1.e-16), xytext=(1.0407e9, 1.e-16),
    arrowprops=dict(arrowstyle="<->"))

plt.text (1.e9, 1.e-3, r"$10^4 g cm^{-3} ")
plt.text (1.e9, 1.e-1, "rho5 = 1")
plt.text (1.e9, 1.e1, "rho5 = 10")

xmin, xmax = plt.xlim ()
ymin, ymax = plt.ylim ()

#plt.ylim (1.e-22, 1.e7)
#plt.xlim (1.e5, 1.e9)

plt.savefig ("enucratio_mod.png", dpi=1000)
