import math
import scipy.misc
import scipy.special

# Calculate exact turbulent enhancmeent factor, as a function of both
# the RMS temperature and density fluctuations on the integral scale L,
# as well as the correlation coefficient between density and temperature.

# Revision history :
# 051919 - Included exact two-body enhancement including non-zero
# density-temperature correlations.

def enhance (n, dT, drho, r):

  sum = 0.

# Following assumes we are calculating the mean turbulent enhancement on the
#  scale r = L. An additional factor must be included when r < L.

  for i in range (0, n):
    if (i % 2 == 0) :
      coeff = math.factorial (n) / (math.factorial (n - i) * math.factorial (i) ) * scipy.special.factorial2 (i - 1) # * (r / L)**(i/3.)
    else:
      coeff = r * drho * math.factorial (n) / (math.factorial (n - i) * math.factorial (i) ) * scipy.special.factorial2 (i) # * (r / L)**((i + 1) /3.)
    term = coeff * dT**i
    print ("i = ", i, " coeff = ", coeff, " term = ", term)
    sum = sum + term
  return sum

# Two-body interactions

# RMS temperature on scale L divided by mean temperature delta T / T_0 
deltaT = 0.08304253333824023

# RMS density on scale L divided by mean density delta rho / rho_0
deltarho = 0.11883278441365815

# Correlation coefficient between density and temperature
rcorr = 0.49037211933846486

#n_T = 4   # p = p
n_T = 23  # C - C
#n_T = 40. # triple-alpha

T9 = 1.
n_T = int (28.05 * T9**(-1./3.) - 2./3.)
print ( "n_T = ", n_T )

#print (enhance (41, x) )
print (enhance (n_T, deltaT, deltarho, rcorr) )
print (enhance (n_T, 0.4, 0., 0.) ) 
print (enhance (n_T, 0.0008, 0., 0.) - 1. )
#print (enhance (40, 0.08) )

