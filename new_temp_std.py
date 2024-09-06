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

 
#GC temperature distribution
  r = temparr/rhoarr
  a = [r[5709], r[5810]] # r[4000], r[4500], r[5298]]
  temp_std = np.std(a, dtype=np.float64)
  temp_dist = temp_std/r[5709]
  #a = np.array(temparr[3490]/rhoarr[3490], temparr[3590]/rhoarr[3590])
  #temp_dist = np.std(a, dtype=np.float64)
  print('Temp distribution =', temp_dist)

