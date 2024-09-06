#Rcrit - strong enhancement regime!!
# -*- coding: utf-8 -*-

#____last update on Jan 12, 2020
#____part of the code baswed on Bob Fisher code and my matlab code
#----@author: Yossef Zenati

#################################################
import scipy.integrate as integrate
import scipy.special as special 
import os
import sys
import csv
import numpy as np 
import matplotlib.pyplot as plt
from decimal import Decimal
import math
#import astropy
#from astropy.constants import G, M_sun


pi = math.pi
enhance = []; xaxis = []; variance = []
enhance2 = []; variance2 = []
imax = 100 #500

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


def turb_enhance (y, temp_pow = 23.):


  mean = integrate.quad(lambda x: (2. * pi)**(-0.5) * (1. + y * x)**temp_pow * math.exp (-x**2./2), -np.inf, np.inf)

  return (mean [0])


L = 80.e5 # 50 ,60 km

#rho9, t8 = 1.e0, 20.0
rho9, t8 = 4.1e-2, 9.2 # He detonation
#cs =  4352.9 * 1.e5 # rho9 = 1.e-2, t8 = 20
mu = 2.34 * 1.6726219e-24 #proton mass
kb = 1.38065*1e-16; #Boltzmann erg/K
gamma = 5./3. #special.gamma ( (n + 1.) / 2.)
cs = (gamma * kb * 1e8*t8 / mu)**0.5
T = t8 * 1.e8
rho = rho9 * 1.e9
#n0 = 28.05 * ((t8/10)**(-1/3)) - 2./3.
xc12 = 0.11 # mass fraction of c12
n0 = (3/(xc12)) - 3 #For He detonation, equation 4 Gronow,S ar el 2020
n = round(n0)
cp = cp (rho9, t8) 
###Fisher
epsilonnuc = []
taunuc = []
xaxis = []
ratio = []

print ("rho9 = ", rho9)
print ("t8 = ", t8)
print ("reaction rate exponent n = ", n)
print ("background sound speed (km/s) = ", cs/1.e5)
print ("integral scale L (km) = ", L / 1.e5)
taunuc0 = cp * T / (n * epsilon (rho9, t8) )
print ("taunuc (s) = ", taunuc0)
print ("Sound speed times nuclear burning time (km) = ", taunuc0 * cs / 1.e5)
######
dToT = [0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.0085,0.009,0.0095,0.0097,0.010,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.085,0.090,0.095,0.100,0.150,0.200,0.250,0.300,0.350,0.400,0.450,0.500,0.550,0.600,0.650,0.700]
dToverTa = np.array(dToT)
#print ("dtt", dToverTa)
#print ("dtt", dToverTa[:5])
#FileRcritRoot = open("RcritvsdtT0" + "_T8_" + str(t8) +  "_Rho9_" + str(rho9) +".dat", "w")
FileRcritRoot = open("RcritvsHedtT0" + "_T8_" + str(t8) +  "_Rho9_" + str(rho9) +".dat", "w")
for ij in range(43):

  #v0 = 2200. * 1.e5 #same
  v0 = 2500. * 1.e5 #He detonation
  imax = 1000 #500
  foundroot = False
  dToverT = dToverTa[ij]
  print("dtt", dToverT)
  #Filertau = open("dToverT_" + str(dToverT) + "_T8_" + str(t8) + "_Rho9_" + str(rho9) + "_L5" + ".dat", "w")
  Filertau = open("Hedet_dToverT_" + str(dToverT) + "_T8_" + str(t8) + "_Rho9_" + str(rho9) + "_L5" + ".dat", "w")
  for i in range (imax):
### Fisher (the same)
    power = -3.3 + float (i / imax) * 3.
    x = 10.**(power)  
    xaxis.append (x)

    epsilon_enhance = epsilon (rho9, t8) * turb_enhance (x, n)
    epsilonnuc.append (epsilon_enhance)
    taunuc_value = cp * T / (n * epsilon_enhance)
    taunuc.append (taunuc_value)
    tausound = L / cs
  
    r = (x / (dToverT) )**3. * L
    #r = (cs * cp * T) / (epsilon_enhance * n * (((n*(n-1.)/2.)*(dToverT**2.)*(x/L)**2.) + 1.))
   
    vr = v0 * (r / L)**(1./3.) #Zenati&Fisher2020
    taueddy = r / vr

    Filertau.write(str(round(r,6)) + "\t" + str(round(vr,6)) + "\t" +str(round(taueddy,6)) + "\t" + str(round(taunuc_value,6)) + "\n")
  
    if ( (tausound / taunuc_value > 1) and (foundroot == False) ) :
#  if ( (taueddy / taunuc_value > 1) and (foundroot == False) ) :
      print ("Root found: x = ", x)
      print ("Turb enhance = ", turb_enhance (x, n) )
      print ("rcrit (km) = ", r / 1.e5)
      FileRcritRoot.write(str(round(r/1.e5,6)) + "\t" + str(dToverT) + "\t" + str(round(taueddy / taunuc_value,6)) + "\n")

      foundroot = True
      
 ##Fisher
#  if ( (tausound / taunuc_value > 1) and (foundroot == False) ) :
#    print ("Root found: x = ", x) 
#    foundroot = True

    ratio.append (taueddy / taunuc_value)
  print ("ratio",ratio )
  Filertau.close()
#  ratio.append (tausound / taunuc_value)

  if (foundroot == False) :
    print ("No root found.")
FileRcritRoot.close()
