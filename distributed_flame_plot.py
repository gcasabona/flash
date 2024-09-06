import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np 
import math

# A short script to calculate the dimensionless enhancement in the mean
# nuclear burning rate, as a function of the dimensionless variable
# delta T / delta T_0 (r / L)^(1/3), where delta T_0 is the RMS variance
# in the temperature field on the scale L.

# Define lists to store data.

pi = math.pi
enhance = []; xaxis = []; variance = []
enhance2 = []; variance2 = []
imax = 1000 # number of points sampled
yminplot = -6. # ylimits
ymaxplot = 10. 

unity = [1.] * imax

def turb_enhance (temp_pow):

  # Some declarations
  xaxis = []; enhance = []; variance = []
  imax = 1000

  for i in range (imax) :

    power = -3.4 + float (i / imax) * 3.
#    power = -3.3 + float (i / imax) * 3.
    y = 10.**(power)  # declare dimensionless argument (dT / T) (r / L)**(1/3)


    mean = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**temp_pow * math.exp (-x**2./2.), -np.inf, np.inf)

    secondmoment = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**(2. * temp_pow) * math.exp (-x**2./2.), -np.inf, np.inf)

    xaxis.append (y)

    variance.append ((secondmoment [0] - mean[0]**2)**(1./2.))
    enhance.append ( (mean [0] - 1.) * 2. / (temp_pow * (temp_pow - 1.) ) )

  return xaxis, enhance # note, not returning variance currently

# Plot the turbulent enhancement versus dimensionless variable.

npp = 4.0
nc12 = 23.0
nta = 41.0

x0, y0 = turb_enhance (npp) # H-burning
x1, y1 = turb_enhance (nc12) # C12-C12 burning
x2, y2 = turb_enhance (nta) # triple-alpha burning
#x3, y3 = turb_enhance (45.) # low temp C12-C12 burning
x3, y3 = turb_enhance (8.) # nu-nubar rate

xta = 0.03

# Define linear transition values
x0lin = [2. / (npp  * (npp  - 1.) )**0.5] * imax
x1lin = [2. / (nc12 * (nc12 - 1.) )**0.5] * imax
x2lin = [xta] * imax
yvertarr = np.logspace (yminplot, ymaxplot, imax)
yvert = yvertarr.tolist ()

ax = plt.axes()
#ax.arrow(1.e-2, 1., 0., 100., width = 0.0001, head_length=0.001, head_width=0.001, fc='b', ec='b')

#xta = 2. / (nta  * (nta  - 1.) )**0.5

# arrows and text labels
ax.annotate("", xy=(2.e-3, 1.), xytext=(2.e-3, 1.e2),
    arrowprops=dict(arrowstyle="<-") )
ax.annotate("", xy=(2.e-3, 1.), xytext=(2.e-3, 1.e-2),
    arrowprops=dict(arrowstyle="<-") )
ax.annotate("", xy=(xta, 1.e8), xytext=(0.02, 1.e8),
    arrowprops=dict(arrowstyle="<-") )
ax.annotate("", xy=(xta, 1.e8), xytext=(0.07, 1.e8),
    arrowprops=dict(arrowstyle="<-") )
plt.text (2.1e-3, 5.e0, "Nonlinear", size = 8)
plt.text (2.1e-3, 0.1, "Weak", size = 8)
plt.text (xta * (1. + 0.0625), 10.**8.25, "Non-Universal", size = 8)
plt.text (xta * (1. - 0.5), 10.**8.25, "Universal", size = 8)
plt.text (2.e-1, 1.e6, r"3$\alpha$", size = 8)
plt.text (1.75e-1, 4.e2, "C12+C12", size = 8)
#plt.text (2.e-1, 5.e-2, "p+p", size = 8)
plt.text (2.e-1, 8.e-2, r"$\nu$", size = 8)

plt.axvspan(xta, 10.**ymaxplot, alpha=0.25, color='red')
plt.axvspan(1.e-3, 0.4, ymin = 10.**yminplot, ymax = 3./8., alpha = 0.25, color = 'blue')

# without nu-nubar rate
#plt.loglog (x0, y0, 'k-', x1, y1, 'k--', x2, y2, 'k-.', linewidth = 2)

# with nu rate and p-p rate
# plt.loglog (x0, y0, 'k-', x1, y1, 'k--', x2, y2, 'k-.', x3, y3, 'k-', linewidth = 2)
# with just nu rate
plt.loglog (x1, y1, 'k-.', x2, y2, 'k--', x3, y3, 'k-', linewidth = 2)
plt.loglog (x0, unity, 'k--', linewidth = 1.0)
plt.loglog (x2lin, yvert, 'k--', linewidth = 1.0)
#plt.loglog (x0, unity, 'k--', x1lin, yvert, 'k--', x2lin, yvert, 'k--', linewidth = 0.5)
plt.xlabel ('RMS Temperature Fluctuation $\delta T / T_0 (r / L)^{1/3}$') 
plt.ylabel ('Burning Enhancment 2 {$[\epsilon_r (T_0) / \epsilon (T_0) ]$ - 1} / (n (n - 1) )')
plt.xlim (1.e-3, 0.4)
plt.ylim (10.**yminplot, 10.**ymaxplot)
#plt.show()

# Write the data to a text file, using numpy savetxt

np.savetxt ('enhancepp.dat', np.array ([np.asarray (x0), np.asarray (y0)])  )
np.savetxt ('enhancec12.dat', np.array ([np.asarray (x1), np.asarray (y1)])  )

plt.savefig ("turbulent_enhancement.png", dpi=1000)

