import yt
import numpy as np 
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from yt.visualization.volume_rendering.transfer_function_helper import TransferFunctionHelper
from yt.visualization.volume_rendering.render_source import ColorTransferFunction
from yt.visualization.volume_rendering.api import Scene, VolumeSource

######################
# Transfer Functions #
######################
bounds = [1e6,1e9]
tf_dens = yt.ColorTransferFunction(np.log10(bounds))
tf_dens.add_layers(5, alpha=[0.0,0.1,1.0,0.1,0.1], colormap='PuBu')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

cmap = plt.get_cmap('autumn_r')
new_cmap = truncate_colormap(cmap, 0.2, 0.8)
def linramp(vals, minval, maxval):
    return (vals.max() - 0.8*vals)/(vals.max() - vals.min())
bounds = [3e-3,3e1]
tf_flam = yt.ColorTransferFunction(np.log10(bounds))
tf_flam.map_to_colormap(-1.0, 0.0, colormap=new_cmap, scale_func=linramp)

######################
#      Camera        #
######################
def rotx_mat(angle):
    return np.array([[1,0,0],
                     [0,np.cos(angle),-np.sin(angle)],
                     [0,np.sin(angle),np.cos(angle)]])
def roty_mat(angle):
    return np.array([[np.cos(angle),0,np.sin(angle)],
                     [0,1,0],
                     [-np.sin(angle),0,np.cos(angle)]])
def rotz_mat(angle):
    return np.array([[np.cos(angle),-np.sin(angle),0],
                     [np.sin(angle),np.cos(angle),0],
                     [0,0,1]])
def zoom_out(cam,progress, start, end):
    #progress: 1-100
    npos = start+progress/100.0*(end-start)
    print(npos)
    cam.set_position(npos)
def rotate(cam, progress, tilt=np.pi/4, camdist=0.1):
    viewangle = tilt
    
    zax = np.array([0,0,1])
    rax = np.array([1,0,0])
    
    rotfactor = 0.01
    l = progress
    rotangle = rotfactor*l*2*np.pi
    
    pos = camdist*roty_mat(rotangle)*rotz_mat(viewangle)*zax
    eps = 0.0001

    pos = camdist*np.dot(rotx_mat(viewangle),np.dot(rotz_mat(rotangle),rax))
    pos_dangle = camdist*np.dot(rotx_mat(viewangle),np.dot(rotz_mat(rotangle+eps),rax))
    northvector = np.cross(pos,pos_dangle-pos)/np.sqrt(np.cross(pos,pos_dangle-pos).dot(np.cross(pos,pos_dangle-pos)))
    cam.set_position(pos)
    cam.north_vector = northvector
    cam.switch_orientation()
