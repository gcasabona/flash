#import matplotlib.pyplot as plt
#import scipy.integrate as integrate
from scipy.stats import moment
import numpy as np
import math

# A short script to calculate the dimensionless enhancement in the mean
# nuclear burning rate, as a function of the dimensionless variable
# delta T / delta T_0 (r / L)^(1/3), where delta T_0 is the RMS variance
# in the temperature field on the scale L, using a numerical integration
# scheme.

# Define lists to store data.

def doublefactorial(n):
     if n <= 0:
         return 1
     else:
         return n * doublefactorial(n-2)

#pi = math.pi
size_list = []

for n in list (range (10) ) :
  size_list.append (int (10**n) )

#size_list = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 1000000000]
moment_list = range (21)

#moment_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]

#mean = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**23. * math.exp (-x**2./2.), -np.inf, np.inf)

for size in size_list :
  randarr = np.random.normal (size = size)

  print ("Size = ", size)
#  print (" Mean = ", np.mean (randarr) )
#  print (" Std = ", np.std (randarr) )
  for n in moment_list:
    cntrlmom = moment (randarr, n) 
    exactmom = doublefactorial (n - 1)
    err = (cntrlmom - exactmom) / exactmom 
    print ("   ", n, " error = ", err)
