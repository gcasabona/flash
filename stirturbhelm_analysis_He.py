###!/usr/bin/python

#====================================================================
#
# stirturbhelm_analysis.py - rtf 4/1/2017
# 
# A short script to plot the global quantities output by the 
# StirTurbHelm problem setup. Usage :
#
# python stirturbhelm_analysis <directory1> <directory2> ...
#
# Where directory1, directory2... is a list of directories in which 
# to process. At least one directory is required. Trailing slashes are
# not required in any of the directory arguments; if one is not present
# in input, it is appended in the script.
#
#====================================================================

import string
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import argparse

# Read in directory list from command line using argparse
parser = argparse.ArgumentParser(description='Analyze StirTurbHelm output.')
parser.add_argument('directory', metavar='N', type=str, nargs='+',
                   help='directory to analyze')
#parser.add_argument('--sum', dest='accumulate', action='store_const',
#                   const=sum, default=max,
#                   help='sum the integers (default: find the max)')

args = parser.parse_args()
directorylist = args.directory

# Verbose flag ( = 1 for verbose output, 0 if not)

verbose = 1

# Open data file for reading

for directory in directorylist:

# Add trailing slash if not present already

  directory = os.path.join(directory, '')

  filename = directory + "stirturbhelm_3d.dat"
  f = open (filename, 'r')

# Inititalize data lists

  time = []
  rho = [] ; keng = []; eint = []; cs = []; etot = []; enuc = []
  c12 = [] ; o16 = [] ; he4 = [] ; temp = []
  momx = [] ; momy = [] ; momz = []
  vort = []

# Read data in line by line

  for line in f:
        if line.strip().startswith("#") :
                if (verbose == 1) :
                        print (line),
        else :

# Split the line of data into fields

                lst = line.split()
          
                time.append (float (lst [0]) )
                rho.append  (float (lst [1]) )
                temp.append (float (lst [2]) )
                he4.append  (float (lst [3]) )
                c12.append  (float (lst [4]) )
                o16.append  (float (lst [5]) )
                momx.append (float (lst [6]) )
                momy.append (float (lst [7]) )
                momz.append (float (lst [8]) )
                etot.append (float (lst [9]) )
                keng.append (float (lst [10]) )
                eint.append (float (lst [11]) )
                cs.append   (float (lst [12]) )
                enuc.append (float (lst [13]) )
                vort.append (float (lst [14]) )

# Convert all lists to np arrays

  timearr = np.array (time)
  rhoarr = np.array (rho)
  temparr = np.array (temp) 
  momxarr = np.array (momx)
  momyarr = np.array (momy)
  momzarr = np.array (momz)
  kengarr = np.array (keng)
  eintarr = np.array (eint)
  etotarr = np.array (etot)
  vortarr = np.array (vort)
  csarr   = np.array (cs)
  enucarr = np.array (enuc)
  c12arr  = np.array (c12)
  o16arr  = np.array (o16)
  he4arr  = np.array (he4)

# Hardwire domain size and mean density
  L = 1.e7
  mass = 1.e6 * L**3.

# Compute derived quantities
  times = np.lib.pad (timearr, (0, 1), 'edge')  
  dt = np.diff (times)
  vrmsarr = (2. * kengarr / rhoarr)**(1./2.)
  csarr = csarr / L**3  # normalize sound speed by hard-wired domain volume
#csarr   = (gamma * (gamma - 1.) * eintarr / rhoarr)**(1./2.)
  macharr = vrmsarr / csarr
  momxarr = momxarr / (rhoarr * csarr)
  momyarr = momyarr / (rhoarr * csarr)
  momzarr = momzarr / (rhoarr * csarr)
#turbdisarr = etotarr - kengarr - enucarr * dt   
  turbdisarr = (etotarr - kengarr - eintarr [0]) / timearr - enucarr # Compute rate per unit time

# Generate plots and write them to the _plots subdirectory

# First check to see if the _plots subdirectory exists
  if not os.path.exists (directory + '_plots'):
    os.makedirs (directory + '_plots')

# Mass vs. time
  plt.plot (timearr, rhoarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('Mass (g)')
  plt.savefig (directory + '_plots/' + 'rho.png', format = 'png')
  plt.clf()

# C12 vs. time
  plt.plot (timearr, c12arr / rhoarr)
  plt.xlabel ('Time (s)') 
  plt.ylabel ('Mean C12 Fraction')
  plt.savefig (directory + '_plots/' + 'c12.png', format = 'png')
  plt.clf()

# O16 vs. time
  plt.plot (timearr, o16arr / rhoarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('Mean O16 Fraction')
  plt.savefig (directory + '_plots/' + 'o16.png', format = 'png')
  plt.clf()

# He4 vs. time
  plt.plot (timearr, he4arr / rhoarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('Mean He4 Fraction')
  plt.savefig (directory + '_plots/' + 'he4.png', format = 'png')
  plt.clf()


# Temp vs. time
  plt.plot (timearr, temparr / rhoarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('Mean temperature (K)')
  plt.savefig (directory + '_plots/' + 'temp.png', format = 'png')
  plt.clf ()


# 3D RMS velocity vs. time
  plt.plot (timearr, vrmsarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('3D RMS Velocity (cm/s)')
  plt.savefig (directory + '_plots/' + 'vrms.png', format = 'png')
  plt.clf()

# cs vs. time
  plt.plot (timearr, csarr)
  plt.xlabel ('Time (s)')
  plt.ylabel ('Sound Speed (cm/s)')
  plt.savefig (directory + '_plots/' + 'cs.png', format = 'png')
  plt.clf()

# Mach number vs. time
  plt.xlabel ('Time (s)')
  plt.ylabel ('3D Mach Number')
  plt.plot (timearr, macharr)
  plt.savefig (directory + '_plots/' + 'mach.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Global X-Momentum Error (Mach number)')
  plt.plot (timearr, momxarr)
  plt.savefig (directory + '_plots/' + 'momx.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Global Y-Momentum Error (Mach number)')
  plt.plot (timearr, momyarr)
  plt.savefig (directory + '_plots/' + 'momy.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Global Z-Momentum Error (Mach number)')
  plt.plot (timearr, momzarr)
  plt.savefig (directory + '_plots/' + 'momz.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Global Vorticity (1/s)')
  plt.plot (timearr, vortarr)
  plt.savefig (directory + '_plots/' + 'vort.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Nuclear energy generation (erg/s)')
  plt.semilogy (timearr, enucarr)
  plt.savefig (directory + '_plots/' + 'enuc.png', format = 'png')
  plt.clf()

  plt.xlabel ('Time (s)')
  plt.ylabel ('Turbulent Dissipation and nuclear energy generation (erg/s)')
  plt.plot (timearr, enucarr / turbdisarr)  # Plot ratio of nuclear energy generation rate to turbulent dissipation
#plt.plot (timearr, enucarr)
  plt.savefig (directory + '_plots/' + 'turbdiss.png', format = 'png')
  plt.clf()
