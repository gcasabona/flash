import numpy as np
import matplotlib.pyplot as plt
import yt
import codecs, json

plt.switch_backend('agg')

# In the 'plotfile' specify the name of the files to plot 
#Open 'plotfile' and read file names 
with open("plotfiles.txt","r") as FL:
    plot_files = FL.readlines()

#specify how the graphs will show like line structure, color, text etc  

line_shape = ['k-'] #, 'k--', 'k-.', 'k:']
txt = ['512^3'] #, '256^3', '128^3', '64^3']
line_width = [1.5] #,1.5,1.5, 2]


temp_max = 2.e9 
temp_min = 1.e8
num_bin = 101

temp_div = (temp_max - temp_min)/num_bin
#print 'div = %e' %(temp_div)

temp_bin = np.arange(temp_min, temp_max, temp_div)
#temp_bin = temp_bin / np.mean(temp_bin)


for i in range(0, len(plot_files)):

    # Load the dataset.

    ds = yt.load(plot_files[i].rstrip('\n'))

    dd = ds.r[:,:,:]
    
    temp = dd["temp"][:]
    cell_mass =  dd['cell_mass'][:]
    
    total_cell_mass = np.sum(cell_mass)
    temp_w_mean = np.sum(temp*cell_mass)/total_cell_mass
    print(temp_w_mean)    

    temp_w_mean_squared = np.sum(temp*temp*cell_mass)/total_cell_mass
    print(temp_w_mean_squared)

    stand_dev_temp = (temp_w_mean_squared - (temp_w_mean)**2)**0.5
    print(stand_dev_temp)

    which_bin = np.digitize(temp, temp_bin) 

    cell_mass_bin = np.zeros(len(temp_bin))
    total_cell_mass = np.sum(cell_mass)

    for j in range(len(temp_bin)):
        condition = which_bin == j
        cell_mass_value =  np.extract(condition, cell_mass)
        #print (cell_mass_value)
        cell_mass_bin[j] = np.sum(cell_mass_value)

    temp_div_normalized = temp_div / temp_w_mean
    temp_bin_normalized = temp_bin / temp_w_mean
    cell_mass_prob = cell_mass_bin / (total_cell_mass*temp_div_normalized)

    plt.loglog (temp_bin_normalized, cell_mass_prob, line_shape[i], linewidth = line_width[i], label = txt[i])
    plt.legend ()
    print (np.trapz(cell_mass_prob, temp_bin_normalized))
    
    x_data = temp_bin_normalized.tolist() 
    y_data = cell_mass_prob.tolist()
    file_path = "/work/05351/casabona/stampede2/run_128_test/new_temp_pdf.json" # "/scratch/04733/tg840329/paper_plot/Temp_pdf.json"
    #json.dump(x_data, codecs.open(file_path, 'aw', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)
    #json.dump(y_data, codecs.open(file_path, 'aw', encoding='utf-8'), separators=(',', ':'), sort_keys=True, indent=4)
plt.xlim(2e-1,2)

#plt.ylabel ('Mass (g)')
plt.ylabel ('Probability Density Function')
plt.xlabel ('Temperature/Mean Temperature')
#plt.xlabel ('Density (g/cm^3)')
#plt.xlabel ('Enuc (erg/s/g)')
#plt.title ('Mass Distribution of Enuc at t = ' + str (ds.current_time.in_units('s')) )
#plt.title ('PDF of Density at t = ' + str (ds.current_time.in_units('s')) )
#plt.title ('PDF of Temperature at t = ' + str (ds.current_time.in_units('s')) )
#plt.title ('PDF of Specific Nuclear Energy Rate at t = ' + str (ds.current_time.in_units('s')) )
plt.savefig('Temperature_pdf.png',format='png',dpi=1000)
