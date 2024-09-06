import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.special as special 
import numpy as np 
from decimal import Decimal
import math

# A short script to calculate the dimensionless enhancement in the mean
# nuclear burning rate, as a function of the dimensionless variable
# delta T / delta T_0 (r / L)^(1/3), where delta T_0 is the RMS variance
# in the temperature field on the scale L.

# Define lists to store data.

pi = math.pi
enhance = []; xaxis = []; variance = []
enhance2 = []; variance2 = []
imax = 100

def cp (rho9, t8):

# cp as defined by eqn. 6 in Woosley, Wunsch, and Kuhlen, 2004
# Note this paper defines cp per 10^8 K; here we use plain cgs
 
  value = 9.1e6 + (8.6e5 * t8 / rho9**(1./3.) ) + (3.0e1 * t8**3. / rho9)
  return value

def epsilon (rho9, t8):

# specific nuclear energy generation rate as defined by Woosley, Wunsch,
# and Kuhlen, 2004

  value = 2.8e13 * (t8 / 7.)**23. * (rho9 / 2.)**3.3
  return value

def epsilon_3a (rho, t8):

  x4 = 1.0
  value = 3.e11 * x4**3. * rho**2. * t8**(-3.) * math.exp (-43.5 / t8)
  return value

def turb_enhance (y, temp_pow = 23.):

# y value is the dimensionless argument (dT / T) (r / L)**(1/3)

  mean = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**temp_pow * math.exp (-x**2./2), -np.inf, np.inf)

  return (mean [0])

# Plot the turbulent enhancement versus dimensionless variable.

#x0, y0 = turb_enhance (4.0) # H-burning
#x1, y1 = turb_enhance (23.) # C12-C12 burning
#x2, y2 = turb_enhance (41.) # triple-alpha

#plt.loglog (x0, y0, 'k-', x1, y1, 'k--', x2, y2, 'k-.',linewidth = 2)
#plt.xlabel ('RMS Temperature Fluctuation $\delta T / \delta T_0 (r / L)^{1/3}$') 
#plt.ylabel ('Burning Enhancment $(\epsilon_r (T_0) / \epsilon (T_0) )$ - 1')
#plt.show()

L = 100.e5 # 10 km

#rho9, t8 = 1.e-2, 20.0
#cs =  4352.9 * 1.e5 # rho9 = 1.e-2, t8 = 20

rho9, t8 = 1.e-2, 15.
cs = 4053. * 1.e5 # rho9 = 1.e-2, T8 = 15.

#rho9, t8 = 1.e-3, 10.
#cs = 2953 * 1.e5 # rho9 = 1.e-3, t8 = 10.

#rho9, t8 = 0.1, 15.
#cs = 5648 * 1.e5 # km/s, rho9 = 0.1, t8 = 15.

T = t8 * 1.e8
rho = rho9 * 1.e9
n = 23.  # c12 + c12
#n = 41. # triple alpha

epsilonnuc = []
taunuc = []
xaxis = []
ratio = []

cp = cp (rho9, t8)

print ("rho9 = ", rho9)
print ("t8 = ", t8)
print ("reaction rate exponent n = ", n)
print ("background sound speed = ", cs)
print ("integral scale L (km) = ", L / 1.e5)
taunuc0 = cp * T / (n * epsilon (rho9, t8) )

print ("taunuc = ", Decimal (taunuc0) )
#print ("taunuc = ", '%.2E' & Decimal (taunuc0) )
print ("Sound speed times nuclear burning time (km) = ", '%.2E' % Decimal (taunuc0 * cs / 1.e5))

#lcrit = pi**(0.5) * cs * taunuc0 / special.gamma ( (n + 1.) / 2.) 
#deltaratio = (pi**(0.5) * cs * cp * T / (n * special.gamma ( (n + 1.) / 2) * L * epsilon (rho9, t8) ) )**(1. / n)

#print ("deltaratio = ", deltaratio)
#print ("lcrit (km) = ", lcrit / 1.e5)

#epsilonnuc = epsilon (rho9, t8) * turb_enhance (0.7)

dToverT = 0.2
v0 = 2200. * 1.e5
imax = 1000
foundroot = False

for i in range (imax):

  power = -3.3 + float (i / imax) * 3.
  x = 10.**(power)  # declare dimensionless argument (dT / T) (r / L)**(1/3)
  xaxis.append (x)

  epsilon_enhance = epsilon (rho9, t8) * turb_enhance (x, n)
#  epsilon_enhance = epsilon_3a (rho, t8) * turb_enhance (x, n) 
  epsilonnuc.append (epsilon_enhance)
  taunuc_value = cp * T / (n * epsilon_enhance)
  taunuc.append (taunuc_value)
  tausound = L / cs

  r = (x / (dToverT) )**3. * L
  vr = v0 * (r / L)**(1./3.)
  taueddy = r / vr

  if ( (taueddy / taunuc_value > 1) and (foundroot == False) ) :
    print ("Root found: x = ", x)
    print ("Turb enhance = ", turb_enhance (x, n) )
    print ("rcrit (km) = ", r / 1.e5)
    foundroot = True
 
#  if ( (tausound / taunuc_value > 1) and (foundroot == False) ) :
#    print ("Root found: x = ", x) 
#    foundroot = True

  ratio.append (taueddy / taunuc_value)
#  ratio.append (tausound / taunuc_value)

if (foundroot == False) :
  print ("No root found.")

plt.loglog (xaxis, ratio)
plt.savefig ("critical_length.png", dpi=1000)

