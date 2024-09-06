# estimate runaway temperature
import math
import numpy as np
import matplotlib.pyplot as plt

#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

rho9 = 1.e-3 # 1.e-2
x12 = 0.25 # 0.5
v0 = 2.e8  # 1.e3 km/s
L = 1.e7   # 100 km

#vr = v0 * (r / L)**(1./3.)
#tedd = r / vr

Tcrit = 2.2e9 # Temperature where enuc = eturb
Tmin  = 3.e8  # Minimum temperature on plot
Tmax =  5.5e9 # Maximum temperature on plot
#r = 1.e5   # 1 km
#v = v0 * (r / L)**(1./3.)

# Peak temperatures of C/O mergers
# M. Dan, S. Rosswog, M. Bruggen, P. Podsiadlowski, 2014
# MNRAS, 438, 14 (arXiv:1308.1667)
# 3.229e8 - 1.0407e9
vertline = np.logspace (-20, 0, num = 100)
temp = np.logspace (math.log10 (Tmin), math.log10 (Tmax), num = 100)
t9   = temp / 1.e9

horizontal = temp / temp

ax = plt.axes()
ax.arrow(0, 0, 0.5, 0.5, head_width=0.05, head_length=0.1, fc='k', ec='k')

enuc = 3.65e45 * rho9**(2.5) * x12**2. * np.exp (-65.894 / t9**(1./3.) )
eturb = v0**3. / L

ratio = enuc / eturb

plt.loglog (temp, ratio, 'k', linewidth = 2.5)
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
plt.text (Tcrit, 1.e6, "Nuclear-Dominated")
ax.annotate("", xy=(Tcrit, 1.e5), xytext=(3.e9, 1.e5),
    arrowprops=dict(arrowstyle="<-"))
ax.annotate("", xy=(Tcrit, 1.e-13), xytext=(1.e9, 1.e-13),
    arrowprops=dict(arrowstyle="<-"))
plt.text (Tcrit, 1.e-12, "Turbulent-Dominated", ha = "right")
plt.text (3.5e8, 1.e-20, "WD Merger Peak Temp")
plt.axvspan(Tcrit, Tmax, alpha=0.25, color='red')
ax = plt.axes()
ax.annotate("", xy=(3.22e8, 1.e-21), xytext=(1.0407e9, 1.e-21),
    arrowprops=dict(arrowstyle="<->"))

xmin, xmax = plt.xlim ()
ymin, ymax = plt.ylim ()

plt.ylim (1.e-22, 1.e7)

plt.savefig ("enucratio.png", dpi=1000)
