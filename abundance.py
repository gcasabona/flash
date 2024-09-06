# computes the distribution of the main elements C12, O16, Ni56 and Mg24
# plots the total mass of each element over time
# computes and plots the burning rate (which includes neutrino cooling) over time
# uses the checkpoint files since abundances and nuclear burning are usually not supplied with the plotfiles

import yt
from yt import derived_field
#from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
#import LibCartesian3D_reset
import glob

flashFolder = './' #  '/work/02418/rkashyap/flash4/super3D/run_149/'
basename = 'stirturbhelm_3d_' # 'super3d_' # sets the prefix of the plotfiles and checkpoint files. Usually 'super3d_' for MHD, 'relax_' for purehydro.
useAllCheckpointFiles = True # set True to use glob.glob to select checkpoint files (by default uses all available files)
endcount = 4 # specify if useAllCheckpointFiles == False, gives the number of the last checkpoint file which is used

# ============================ constructs list of filenames if script is run on its own, otherwise supplied by allPlots.py !
def getfilenames(useAllCheckpointFiles,endcount):
   if useAllCheckpointFiles == True:
      chkFilenames_own = glob.glob(flashFolder + basename + 'hdf5_chk_[0-9][0-9][0-9][0-9]')
      chkFilenames_own.sort()
   else:
      chkFilenames_own = []
      for n in range(0,endcount + 1):
         filename = flashFolder + basename + 'hdf5_chk_%04d' % n
         chkFilenames_own.append(filename)
   return chkFilenames_own
 
   print ( 'chkFilenames_own =', chkFilenames_own)
# =========================== beginning of definition of main
@derived_field(name='CellVolume', sampling_type='cell', units = 'cm**3')
def CellVolume(field, data):
    if data['dx'].size == 1:
        try:
            return data['dx'] * data['dy'] * data['dz'] * \
                np.ones(data.ActiveDimensions, dtype='float64')
        except AttributeError:
            return data['dx'] * data['dy'] * data['dz']
    return data["dx"] * data["dy"] * data["dz"]

@derived_field(name='CellMass', sampling_type='cell', units = 'g')
def CellMass(field,data):
        return data['CellVolume']*data['dens']

def main(filenames):
	#print ('Executing abundance.py')
    print ("filenames = ", filenames)
    totalMass = np.zeros(len(filenames))
    he4Mass = np.zeros(len(filenames))
    c12Mass = np.zeros(len(filenames))
    o16Mass = np.zeros(len(filenames))
    ne20Mass = np.zeros(len(filenames))
    mg24Mass = np.zeros(len(filenames))
    si28Mass = np.zeros(len(filenames))
    s32Mass = np.zeros(len(filenames))
    ar36Mass = np.zeros(len(filenames))
    s32Mass = np.zeros(len(filenames))
    ca40Mass = np.zeros(len(filenames))
    ti44Mass = np.zeros(len(filenames))
    cr48Mass = np.zeros(len(filenames))
    ni56Mass = np.zeros(len(filenames))
    fe52Mass = np.zeros(len(filenames))

    time = np.zeros(len(filenames))
    burningRate = np.zeros(len(filenames))

#    conversionfactor = 1.98855e33 # converts from gram to solar mass

    for n in range(len(filenames)):
       pf = yt.load(filenames[n])
       time[n] = pf.current_time
       d = pf.all_data()
       totalMass[n] = (d['CellMass']).sum()
       he4Mass[n] = (d['CellMass'] * d['he4 ']).sum()
       o16Mass[n] = (d['CellMass'] * d['o16 ']).sum()
       c12Mass[n] = (d['CellMass'] * d['c12 ']).sum()
       ne20Mass[n] = (d['CellMass'] * d['ne20']).sum()
       mg24Mass[n] = (d['CellMass'] * d['mg24']).sum()
       si28Mass[n] = (d['CellMass'] * d['si28']).sum()
       s32Mass[n] = (d['CellMass'] * d['s32 ']).sum()
       ar36Mass[n] = (d['CellMass'] * d['ar36']).sum()
       ca40Mass[n] = (d['CellMass'] * d['ca40']).sum()
       ti44Mass[n] = (d['CellMass'] * d['ti44']).sum()
       cr48Mass[n] = (d['CellMass'] * d['cr48']).sum()
       ni56Mass[n] = (d['CellMass'] * d['ni56']).sum()
       fe52Mass[n] = (d['CellMass'] * d['fe52']).sum()

       burningRate[n] = (d['CellVolume'] * d['enuc']).sum()

    conversionfactor = totalMass [0]
    he4Mass /= conversionfactor
    o16Mass /= conversionfactor
    c12Mass /= conversionfactor
    ne20Mass /= conversionfactor
    mg24Mass /= conversionfactor
    si28Mass /= conversionfactor
    s32Mass /= conversionfactor
    ar36Mass /= conversionfactor	
    ca40Mass /= conversionfactor
    ti44Mass /= conversionfactor
    cr48Mass /= conversionfactor
    ni56Mass /= conversionfactor
    fe52Mass /= conversionfactor

# Output text diagnostics 

    print ("Time (s) = ", time)
    print ("He4 abundance = ", he4Mass)
    print ("C12 abundance = ", c12Mass)
    print ("O16 abundance = ", o16Mass)
    print ("Ne20 abundance = ", ne20Mass)
    print ("Mg24 abundance = ", mg24Mass)
    print ("Si28 abundance = ", si28Mass)
    print ("S32  abundance = ", s32Mass)
    print ("Ar36 abundance = ", ar36Mass)
    print ("Ca40 abundance = ", ca40Mass)
    print ("Ti44 abundance = ", ti44Mass)
    print ("Cr48 abundance = ", cr48Mass)
    print ("Fe52 abundance = ", fe52Mass)
    print ("Ni56 abundance = ", ni56Mass)

# Normalization check

    checknorm = he4Mass [0] + c12Mass [0] + o16Mass [0] + ne20Mass [0] +\
                 mg24Mass [0] + si28Mass [0] + s32Mass [0] + ar36Mass [0] +\
                 ca40Mass [0] + ti44Mass [0] + cr48Mass [0] + fe52Mass [0] +\
                 ni56Mass [0]

# Generate plots

    plt.figure(1)
    #plt.subplot(221)
    plt.plot(time,o16Mass,label='O16')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)
	#plt.subplot(222)
    plt.plot(time,c12Mass,label='C12')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)
	

	#plt.subplot(22)
    plt.plot(time, ne20Mass,label='Ne20')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

#plt.subplot(222)
    plt.plot(time, si28Mass,label='Si28')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

	#plt.subplot(222)
    plt.plot(time, ar36Mass,label='Ar36')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

	#plt.subplot(223)
    plt.plot(time,mg24Mass,label='Mg24')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

	#plt.subplot(222)
    plt.plot(time, ca40Mass,label='Ca40')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

	#plt.subplot(224)
    plt.plot(time,ni56Mass,label='Ni56')
    plt.xlabel('time [s]')
    plt.ylabel('mass [m_solar]')
    plt.legend(loc=0)

	#plt.plot(time,n14Mass,label='N14')
	#plt.legend(loc=4)

    plt.plot(time,ne20Mass,label='Ne20')
    plt.legend(loc=0)

    plt.plot(time,fe52Mass,label='Fe52')
    plt.legend(loc=0)

    plt.tight_layout()
    plt.savefig('abundance_elements.png')

    plt.figure(2)
    plt.plot(time,burningRate,label='Burning Rate')
    plt.xlabel('time [s]')
    plt.ylabel('burning rate [erg/s]')

    plt.tight_layout()
    plt.savefig('abundance_burning.png')

#========================================== end of definition of main
if __name__ == '__main__':
   chkFilenames_own = getfilenames(useAllCheckpointFiles,endcount)
   print ("chkFilenames_own dunder = ", chkFilenames_own)
   main(chkFilenames_own)
