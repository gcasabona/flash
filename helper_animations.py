import yt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
from matplotlib import rc_context
import matplotlib.animation as animation
import tqdm
import pickle
import numpy as np
import time
import time
from .helper_plots import slice_plot,slice_create,PlotObj
import matplotlib.pyplot as plt
import shutil
import fileinput
import numpy

def zoomfactor(i,zoomdct=None):
    if zoomdct is None:
        return 1.0
    minz = min(zoomdct["start"][0],zoomdct["end"][0])
    maxz = max(zoomdct["start"][0],zoomdct["end"][0])
    z = (zoomdct["end"][0]-zoomdct["start"][0])/(zoomdct["end"][1]-zoomdct["start"][1])*(i-zoomdct["start"][1])+zoomdct["start"][0]
    z = max(minz,min(maxz,z))
    return z


def create_animation_timeseries(path,field,save=False,mode="mp4",width=2500,zlim=None,cmap="viridis",log=True,interval=100,zoffset=None,debug=False,
                                zoom=None,imagetype="slice",projmethod="integrate",Nframes=None):
    """zoom: Time dependent zoom for slice IDs. E.g.: zoom={"start":[1,0],"end":[5,50]} will zoom from factor 1 to 5 from frame 0 to frame 50."""
    ts = yt.load(path)


    if type(field)!=tuple: # assume field to be a flash/custom flash derived field.
        field = ("flash",field)

    comm = yt.communication_system.communicators[-1]
    ranks = comm.size
    rank = comm.rank

    if Nframes is None:
        Nframes = len(ts)

    if rank==0:
        print("%i frames will be rendered."%Nframes)

    if ranks==1:
        pbar = tqdm.tqdm(total=Nframes) # terminal version
    def compute_slice(i,width=2500):
        ds = ts[i]
        fparr = "tmp/"+save+"_frame"+str(i).zfill(6)+".hdf5"

        pltobj = PlotObj()
        width = width*zoomfactor(i,zoom)
        pltobj.image_create(ds,field,width=(width,"km"),imagetype=imagetype,method=projmethod,debug=debug)
        pltobj.image_create(ds,("flash","flam"),width=(width,"km"),imagetype=imagetype,method=projmethod,debug=debug,fieldid=1)
        pltobj.image_create(ds,("flash","density"),width=(width,"km"),imagetype=imagetype,method=projmethod,debug=debug,fieldid=2)
        pltobj.save(fparr)
        if ranks==1:
            pbar.update(1)
        else:
            print("Frame %i computed."%i)


    comm.barrier()
    t1 = time.time()
    for t in ts.piter(dynamic=False):
        i = int(t.basename.split("_")[-1])
        if i>=Nframes:
            break
        compute_slice(i,width=width)
    t2 = time.time()

    comm.barrier()

    def animate(i):
        fparr = "tmp/"+save+"_frame"+str(i).zfill(6)+".hdf5"
        pltobj = PlotObj()
        pltobj.load(fparr)
        data = pltobj.imarrs["0"]
        ax.clear()
        im = pltobj.plot_ax(ax,timestamp=timestamp,maskfield=maskfield,backgroundfield=backgroundfield,flamcontour=True)
        print("Rendered frame %i."%i)
        return im

    timestamp = True
    maskfield = ("flam","above",0.1)
    backgroundfield = ("flash","density")


    #if yt.is_root(): ;; does not work
    #print("mpirank_debug",rank,yt.is_root())
    if rank==0:
        print("Computation of slices took %.1f seconds."%(t2-t1))

        fparr = "tmp/"+save+"_frame"+str(i).zfill(6)+".hdf5"
        pltobj = PlotObj()
        pltobj.load(fparr)
        fig,ax = pltobj.plot(timestamp=timestamp,maskfield=maskfield,backgroundfield=backgroundfield,flamcontour=True)

        ani = FuncAnimation(fig, animate, frames=Nframes,interval=interval)
        if save:
            outputwritten = False
            with rc_context({'mathtext.fontset': 'stix'}):
                if "mp4" in mode:
                    writer = animation.writers['ffmpeg'](fps=1000/interval,bitrate=-1,codec="libx264")
                    ani.save('./results/'+save+'.mp4',dpi=300,writer=writer)
                    outputwritten = True
                if "webm" in mode:
                    writer = animation.writers['ffmpeg'](fps=1000/interval,bitrate=-1,codec="libvpx-vp9",
                                                         extra_args=["-b:v","2M"])
                    ani.save('./results/'+save+'.webm',dpi=300,writer=writer)
                    outputwritten = True
                if "html" in mode:    
                    writer = animation.writers['html']()
                    fp = "./results/"+save+'.html'
                    ani.save(fp,writer=writer)
                    ###
                    # bugfix for matplotlib.animation behavior for subdired filepath. 
                    with fileinput.FileInput(fp, inplace=True, backup='.bak') as file:
                        for line in file:
                            print(line.replace("./results/", "./"), end='')
                    ###
                    outputwritten = True
                if not outputwritten:
                    print("Unknown 'mode' for output.")
        return ani

    if ranks==1:
        pbar.close()   
    return None # return None if not root process

def create_animation_timeseries_mpl(path,field,save=False,mode="movie",width=2500,zlim=None,cmap="viridis",log=True,interval=100):
    """ Animation simplified with matplotlib.animation. However not parallizable because of it. Thus depreciated.
        save: save standalone animation to disk (default: False).
        modes: 'movie' and 'html'
        path: Path of plot files matching wildcards, e.g. 'gcd_hdf5_plt_cnt_*'"""
    ts = yt.load(path)

    zoffset = -50
    plot = yt.SlicePlot(ts[0], 'x', field,width=(width,"km"))
    plot.set_center((0,zoffset+width/2),"km")
    if zlim:
        plot.set_zlim(field, *zlim)
    plot.set_log(field,True)
    plot.set_cmap(field,cmap)
    fig = plot.plots[field].figure

    Nframes = len(ts)
    
    pbar = tqdm.tqdm(total=Nframes) # terminal version
    #pbar = tqdm.notebook.tqdm(total=Nframes) # modern version
    def animate(i):
        ds = ts[i]
        plot._switch_ds(ds)
        pbar.update(1)

    ani = FuncAnimation(fig, animate, frames=Nframes,interval=interval)
    

    if save:
        with rc_context({'mathtext.fontset': 'stix'}):
            if mode=="movie":
                writer = animation.writers['ffmpeg'](fps=30,bitrate=-1,codec="libx264")
                ani.save('./results/'+save+'.mp4',dpi=150,writer=writer)
            else:    
                writer = animation.writers['html']()
                ani.save('./results/'+save+'.html',writer=writer)
    pbar.close()   
    return ani
