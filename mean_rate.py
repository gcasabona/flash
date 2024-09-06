import matplotlib.pyplot as plt
import scipy.integrate as integrate
import numpy as np
import math

# A short script to calculate the dimensionless enhancement in the mean
# nuclear burning rate, as a function of the dimensionless variable
# delta T / delta T_0 (r / L)^(1/3), where delta T_0 is the RMS variance
# in the temperature field on the scale L, using a numerical integration
# scheme.

# Define lists to store data.

pi = math.pi
enhance = []; xaxis = []; variance = []

y = 0.007714211121331792 # delta T / T_0 (r / L)^(1/3)

mean = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**23. * math.exp (-x**2./2.), -np.inf, np.inf)

print (mean [0] - 1.)
