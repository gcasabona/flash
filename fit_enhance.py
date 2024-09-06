#!/usr/bin/env python

# Fit the turbulent nuclear enhancement rate to a simplified functional
#  form. Uses the LMFIT non-linear least-squares minimiation and curve
#  fitting package for Python (https://lmfit.github.io/lmfit-py/).

# RTF. Last modified 6/6/2019.
 
import matplotlib.pyplot as plt
from numpy import loadtxt
import numpy as np
import math

from lmfit import Model, Parameters

# Load the computed turbulent enhancement factor from a file.

data = loadtxt('enhancec12.dat')
print (data)

#n = 4. # pp
n = 23. # c12
econst = math.e

def custom(x, A, x0, alpha, beta):
    """Define a broken power-law model."""
    return A * x**alpha * (1. + (x/x0))**beta 

def logcustom(x, A, x0, alpha, beta, C):
    """Define a log broken power-law model."""
    return np.log (A * x**alpha * (1. + econst**(beta * (x - x0)) )  ) + C

x = data [0, :]
y = data [1, :]

logy = np.log (y)
logx = np.log (x)

cmodel = Model(logcustom)
params = Parameters ()

#define parameters of the fit with constraints
params.add ('A', value = 0.05, vary = True, min = 0, max = 1.)
params.add ('x0', value = 0.0625, vary = True, min = 0.0, max = 1.)
params.add ('alpha', value = 2., vary = False, min = 0.0, max = 10.)
params.add ('beta', value = 30, vary = True, min = 1., max = 50.)
params.add ('C', value = 5.5, vary = True, min = 0., max = 10.)

#result = cmodel.fit(logy, x=x, A = 1., x0 = 0.05, alpha  = 2.,beta = n / 3., nan_policy='propagate')
result = cmodel.fit(logy, params, x=x)

print(result.fit_report())

diff = logy - result.best_fit
fracdiff = abs (diff) / logy

print ("Max fractional diff = ", max (fracdiff) ) 

plt.plot(x, logy, 'k--')
#plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()
# <end examples/doc_model_gaussian.py>
