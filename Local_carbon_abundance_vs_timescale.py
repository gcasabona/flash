import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
    
file_name = ['tburn.dat']
time = []
t_burn = []
c12_abundance = []

#GC

f = open (file_name[0], 'r')
counter = 0
for line in f:
  counter += 1
  if (counter> 0 and counter <2527):
    if (counter%2):
      #print(line)
      lst = line.split()
      #print(lst)
      time.append (float (lst [0]))
      c12_abundance.append(float (lst [2]))
    
    else:
      print(line)
      lst = line.split()
      print(lst)
      t_burn.append(float(lst [0]))

      print(time)
      print(min(t_burn))
      print(max(t_burn))
      plt.semilogy(time,t_burn)
      plt.savefig('tburn_vs_time.png',format ='png',dpi = 300)
      plt.show()
      print(time[0])
      new_time = np.ones(len(time))*time[0]
      time = time - new_time
    

      eddy_length_scale = 1.0  # unit km
      L_box = 100.0            # unit km
      global_rms_velocity =  2000.0 #unit km/s
      eddy_time = eddy_length_scale/(global_rms_velocity*(eddy_length_scale/L_box)**(1.0/3.0))  # r/(v_o*(r/L)**1/3)
      eddy_time = np.ones(len(t_burn))*eddy_time
      np.mean(eddy_time)
      fig, axis_1 = plt.subplots()
      axis_2 = axis_1.twinx()
      color = 'tab:blue'
      axis_1.set_xlabel('time (s)', fontsize = 12)
      axis_1.set_ylabel('t_eddy/t_burn', color = color)
      axis_1.set_ylabel(r'$t_{edd}$ / $t_{burn}$', fontsize = 12 ,color = color)
      axis_1.semilogy(time,eddy_time/np.asarray(t_burn),color=color)
      axis_1.tick_params(axis='y',labelcolor = color)
      axis_1.legend(loc = 'center left')
      color = 'tab:red'
      axis_2.set_ylabel('c12_abundance', color = color)
      axis_2.set_ylabel(r'min ($X_{12}$ abundance)',fontsize =12, color = color)
      axis_2.plot(time, np.asarray(c12_abundance),color = color)
      axis_2.tick_params(axis='y',labelcolor = color)
      axis_2.legend(loc = 'center left')
      fig.savefig('local_timescale_vs_local_c12.png',format = 'png', dpi = 1000)
      plt.show()
