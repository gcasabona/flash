###########################################################
## 
## analysis_pdf.py, rtf 4/14/17
##
## A short yt-based script to generate a one-point PDF of
## simulation data, and to fit this to a prescribed 
## function.
##
#############################################################

import yt
from yt.units import gram, second
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
from scipy.optimize import leastsq
from lmfit.models import PowerLawModel, ExpressionModel, Model
import argparse
import glob

# Define a power-law fitting function, with y-offset.
#powerlaw = lambda x, amp, off, coeff, index : amp * ((x - off)**index) + coeff
#powerlaw = lambda x, amp, coeff, index : amp * x**index + coeff
#fitfunc = lambda p, x: p[0]*(x**p[1]) + p [2] 
#fitfunc = lambda p, x: p[0] * (( x - p[1] )**p[2]) + p[3]
#errfunc = lambda p, x, y: (y - fitfunc(p, x) )
def mypowerlaw (x, amp, index, xoff, yoff) :
  return amp * (x - xoff)**index + yoff

# Read in directory list from command line using argparse
parser = argparse.ArgumentParser(description='Analyze StirTurbHelm output.')
parser.add_argument('filelist', metavar='N', type=str, nargs='+',
                   help='list of files to analyze')
parser.add_argument('field', metavar='Field to analyze PDF', type=str, nargs=1,
                   help='field to obtain PDF')
args = parser.parse_args()
filelist = args.filelist
field = args.field

# Load the dataset.

for file in filelist:

#  print "Searching files in ", directory + "stirturbhelm*hdf5*" 
#  filelist = glob.glob (directory + "stirturbhelm*hdf5*")
#  print "length of filelist = ", len (filelist)
#  print "filelist = ", filelist
#  filelist = filelist.sort()
#  print "sorted filelist = ", filelist
#  print filelist [-1]
#  lastfile = filelist [-1]
  ds = yt.load(file)

# Create a data object that represents the whole box.

  ad = ds.all_data()

# Generate the PDF, weighted by cell mass. This produces density with individual
# bins of mass, which we later normalize to the total mass on the domain.

#plot = yt.ProfilePlot(ad, "temperature", ["cell_mass"], weight_field=None)
#plot = yt.ProfilePlot(ad, "density", ["cell_mass"], weight_field=None)
  plot = yt.ProfilePlot(ad, "temp", ["cell_mass"], weight_field=None)

# Define region encompassing the total domain
  domain = ds.region(ds.domain_center, ds.domain_left_edge, ds.domain_right_edge)

# Get the total volume and density on the domain for use in normalizing the
# PDF
  volume = domain.sum ("cell_volume")
  mean_density = domain.mean ("density")

# Pull in the object corresponding to the PDF plot.

  profile = plot.profiles [0]

# Extract the raw data for the PDF plot.
  yt_xdata = profile.x
  yt_ydata = profile['cell_mass']

  xdata = yt_xdata.v
  ydata = yt_ydata.v / (mean_density * volume) # normalized to total mass

# Test data
#x = [1., 2., 3., 4., 5., 6.]
#y = [0.8, 4.8, 9.1, 15.3, 25.9, 38.]
#xdata = np.asarray (x)
#ydata = np.asarray (y)

  xdata = xdata [np.where(np.isfinite (xdata ) ) ]
  ydata = ydata [np.where(np.isfinite (ydata ) ) ]

  max_x = np.amax (xdata)
  max_y = np.amax (ydata)

  xdata = xdata / max_x
  ydata = ydata / max_y 

  print ("xdata = ", xdata)
  print ("ydata = ", ydata)
# Fit the PDF to the powerlaw form specified above.

#model = Model (mypowerlaw)
  #GC model = GaussianModel ()
  #GC params = model.guess (ydata, x=xdata)
  #GC result = model.fit (ydata, pars, x = xdata)
#result = model.fit (ydata, x=xdata, exponent = -1., amplitude = 1.)
#result = model.fit (ydata, x=xdata)

#GC temp std
  d = np.array([xdata, ydata])
  temp_std = np.std(xdata)
  print("temp distribution = ", temp_std)
'''
#model = ExpressionModel('amp * (x - xoffset)**index + yoffset')
#result = model.fit (ydata, x = xdata, amp = 1., index = 1., xoff = 0.5, yoff= 0.5)
#print (result.fit_report() )
#print (result.eval_uncertainty () )

#model.plot_fit()

#  plt.loglog(xdata * max_x, ydata * max_y,         'bo')
  plt.loglog(xdata * max_x, result.init_fit * max_y, 'k--')
#  plt.loglog(xdata * max_x, result.best_fit * max_y, 'r-', linewidth = 2)
#  plt.show()
#  plt.clf()

# Plot the data along with the best-fit.

#plt.loglog (max_x * xdata, max_y * ydata, linewidth = 2)
#plt.loglog(max_x * xdata, max_y * powerlaw(xdata, amp, off, coeff, index), 'g--', label='fit')
#plt.loglog(max_x * xdata, max_y * powerlaw(xdata, amp, coeff, index), 'g--', label='fit')
  plt.loglog (xdata * max_x, ydata * max_y, linewidth = 2)
#plt.ylabel ('Mass (g)')
  plt.ylabel ('Probability')
  plt.xlabel ('Temperature (K)')
#plt.xlabel ('Density (g/cm^3)')
#  plt.xlabel ('Enuc (erg/s/g)')
#plt.title ('Mass Distribution of Enuc at t = ' + str (ds.current_time.in_units('s')) )
#plt.title ('PDF of Density at t = ' + str (ds.current_time.in_units('s')) )
  plt.title ('PDF of Temperature at t = ' + str (ds.current_time.in_units('s')) )
#  plt.title ('PDF of Specific Nuclear Energy Rate at t = ' + str (ds.current_time.in_units('s')) )
  plt.legend ()
#  plt.show ()

# Save the PDF.
# Optionally, give a string as an argument
# to name files with a keyword.
  plt.savefig('driventurb_pdf_' + file[-4:] )

'''
