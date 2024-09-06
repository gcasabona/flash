import yt
import numpy as np
# Load the dataset (replace 'data_filename' with your file's name)
ds = yt.load('stirturbhelm_3d_hdf5_chk_0083')
# Access the data for the field "temp"
ad = ds.all_data()
temp_data = ad["temp"]
# Compute the standard deviation of the "temp" field
std_temp = np.std(temp_data)
# Compute the mean value of the "temp" field
mean_temp = np.mean(temp_data)
print(f"RMS value of 'temp': {std_temp}")
print("RMS / mean temp", std_temp / mean_temp)
