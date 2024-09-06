import yt
import matplotlib.pyplot as plt
import numpy as np
import glob

# Parameters for analysis; names of files and colors of curves
#filenames = ['run64', 'run128', 'run256', 'run512']

def periodic_distance (a, b, L):

  diff = a - b
  for i in range (len (diff) ):
    if (abs (diff [i]) > L / 2.) :

  distance = np.linalg.norm (a - b, 2)  
  return distance
 
# === Parameters utilized in execution ===

filenames = glob.glob ("_run256/*0")
filenames.sort() 
print (filenames)

# data field to analyze
#data_field = 'temperature'
data_field = 'pressure'
#data_field = 'enuc'
#data_field = "density"

logscales = {
  "radius" : False,
  data_field : False
}

units = {
  "radius" : "km",
  data_field : "K" # temperature
}

nbins = 10

# === end of parameters ===

# Declare empty lists to be used later
datalist = []; xlist = []; ccdflist = []

profiles = []

c_previous = [0, 0, 0]
c_previous_vec = np.asarray (c_previous)

for filename in filenames :
  ds = yt.load(filename)
  v, c = ds.find_max ("enuc")
  print ("Centered at v = ", v)
  print ("Centered at c = ", c)
# Print the distance between this hot spot maximum and the previous one. These should be relatively small if we are tracking the same hot spot
  print ("distance = ", np.linalg.norm (c - c_previous_vec, 2) / 1.e5 )
  my_sphere = ds.sphere(c, (5, "km") )
#  profiles.append(yt.create_profile (my_sphere, ["radius"], fields = [data_field], weight_field = "cell_mass", n_bins=nbins, logs = logscales, units = units) )
  profiles.append(yt.create_profile (my_sphere, ["radius"], fields = [data_field], weight_field = "cell_mass", n_bins=nbins, logs = logscales) )
#  detplot = yt.ProfilePlot (my_sphere, "radius", "temperature")
  c_previous = c

plot = yt.ProfilePlot.from_profiles(profiles)

plot.save()
