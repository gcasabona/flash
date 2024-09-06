import yt
from yt.visualization.plot_window import get_window_parameters
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,Normalize
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.patheffects as PathEffects
import numpy as np

from .helper_ytfields import fieldplotprops,flameprops
fpps = fieldplotprops # get from helpers_ytfields

import h5py

class PlotObj():
    """Holds relevant information for Slice/Projection/VolumeRender plots."""
    def __init__(self,width=(5000,"km"),center=[0.0,0.0,0.0]):
        self.params = {}
        self.params["center"] = center
        if type(width)!=tuple:
            raise ValueError("Only accepting tuple (value,unitstring) for width.")
        self.params["width"] = width
        self.imarrs = {}
    def image_create(self,ds,field,fieldid=0,imagetype="slice",**kwargs):
        kwargs["center"] = kwargs.get("center",self.params["center"])
        kwargs["width"] = kwargs.get("width",self.params["width"])
        keyid = ""
        if fieldid>0:
            keyid += "_%i"%fieldid

        if imagetype=="slice":
            self.params["type"+keyid] = "SlicePlot"
        elif imagetype=="projection":
            self.params["type"+keyid] = "ProjectionPlot"
        else:
            raise ValueError("Not implemented imagetype.")

        imarr = self.frb_create(ds,field,imagetype=imagetype,**kwargs)

        self.imarrs["%i"%fieldid] = imarr
        self.params["field"+keyid] = field
        self.params["field_display_name"+keyid] = self.params["field"][1]
        try:
            displayname = ds.field_info[field].display_name
            if displayname:
                self.params["field_display_name"+keyid] = displayname
        except AttributeError: # fall back onto field name
            pass # just accept the priorly set fall-back value.
        self.params["units"+keyid] = ds.field_info[field].units # units of plotted quantity
        self.params["snapshot_time"] = ds.current_time

    def frb_create(self,ds,field,debug=False,axis=None,fieldid=0,resolution=1600,imagetype="slice",
                   method="integrate",**kwargs):
        """ Creates a slice and returns 2D numpy array of such.
            fieldid: Specify as integer >=1 to save as secondary field under given identifier."""
        if axis is None: # TODO: Adjust for 1D/2D
            axis = yt.funcs.fix_axis("x",ds) # default for 3D
        if debug:
            resolution = 100 # lower resolution for faster debugging
        (bounds, center, display_center) = get_window_parameters(axis, kwargs["center"], kwargs["width"], ds)

        keyid = ""
        if fieldid>0:
            keyid += "_%i"%fieldid

        self.params["extent"] = bounds
        data_source = ds.all_data()
        data_source.set_field_parameter("flameprops",flameprops)
        if imagetype=="slice":
            p = ds.slice(axis, center[axis],center=center,data_source=data_source)
        else:
            p = ds.proj(field,axis,center=center,data_source=data_source,method=method)
        frb = p.to_frb(kwargs["width"], resolution)
        return frb[field]



    def proj_create(self,ds,field,fieldid=0,**kwargs):
        data_source = ds.all_data()
        data_source.set_field_parameter("flameprops",flameprops)
        p = yt.ProjectionPlot(ds, "x", field, width=kwargs["width"],data_source=data_source)
        return p.frb[field]
    #def slice_create(self,ds,field,debug=False,axis=None,fieldid=0,**kwargs):
    #    """ Creates a slice and returns 2D numpy array of such.
    #        fieldid: Specify as integer >=1 to save as secondary field under given identifier."""
    #    kwargs["center"] = kwargs.get("center",self.params["center"])
    #    kwargs["width"] = kwargs.get("width",self.params["width"])
    #    #kwargs["extent"] = [kwargs.get("width",self.params["width"])

    #    resolution = 1600
    #    if axis is None: # TODO: Adjust for 1D/2D
    #        axis = yt.funcs.fix_axis("x",ds) # default for 3D
    #    if debug:
    #        resolution = 100 # lower resolution for faster debugging
    #    (bounds, center, display_center) = get_window_parameters(axis, kwargs["center"], kwargs["width"], ds)

    #    keyid = ""
    #    if fieldid>0:
    #        keyid += "_%i"%fieldid


    #    self.params["extent"] = bounds
    #    data_source = ds.all_data()
    #    data_source.set_field_parameter("flameprops",flameprops)
    #    slc = ds.slice(axis, center[axis],center=center,data_source=data_source)
    #    frb = slc.to_frb(kwargs["width"], resolution)
    #    self.imarrs["%i"%fieldid] = frb[field]
    #    self.params["type"] = "SlicePlot"
    #    self.params["field"+keyid] = field
    #    self.params["field_display_name"+keyid] = self.params["field"][1]
    #    try:
    #        displayname = ds.field_info[field].display_name
    #        if displayname:
    #            self.params["field_display_name"+keyid] = displayname
    #    except AttributeError: # fall back onto field name
    #        pass # just accept the priorly set fall-back value.
    #    self.params["units"+keyid] = ds.field_info[field].units # units of plotted quantity
    #    self.params["snapshot_time"] = ds.current_time

    def plot_ax(self,ax,timestamp=False,scalebar=500,maskfield=None,backgroundfield=None,flamcontour=False,fpp={},bfpp={},info_height=0.15,
                scalebar_pos=None):
        if "type" not in self.params:
            raise KeyError("No object to plot.")
        field = self.params["field"]
        if not fpp:
            fpp = fieldplotprops.get(field,fpp)
        print("fpp",fpp)
        imarr = self.imarrs["0"]

        ims = {}

        if maskfield is not None:
            mid,imarr_maskfield = self.get_imarr(maskfield[0])
            if maskfield[1]=="above":
                imarr = np.ma.masked_where(imarr_maskfield<maskfield[2],imarr)
            else:
                imarr = np.ma.masked_where(imarr_maskfield>maskfield[2],imarr)
        if backgroundfield is not None:
            mid,imarr_bgf = self.get_imarr(backgroundfield)
            fpp_bgf = fieldplotprops.get(backgroundfield,{})
            if bfpp:
                fpp_bgf = bfpp
            norm = get_norm(fpp_bgf)
            cmap = get_cmap(fpp_bgf)
            ims["bgf"] = {"im": ax.imshow(imarr_bgf,norm=norm,origin="bottom",extent=self.params["extent"],cmap=cmap),
                          "field":backgroundfield,"field_display_name":self.params.get("field_display_name_%i"%mid,""),
                          "units":self.params.get("units_%i"%mid,"")}
        if flamcontour:
            self.add_flamcontour(ax) 
        norm = get_norm(fpp)
        cmap = get_cmap(fpp)
        ims["main"] = {"im": ax.imshow(imarr,norm=norm,origin="bottom",extent=self.params["extent"],cmap=cmap),
                       "field":field,"field_display_name":self.params.get("field_display_name",""),
                       "units":self.params.get("units","")}
        ax.set_xlabel(r"p$_1$ [%s]"%self.params["width"][1])
        ax.set_ylabel(r"p$_2$ [%s]"%self.params["width"][1])
        if timestamp:
            plot_time(ax,self.params["snapshot_time"],center=[0.85,info_height])
        if scalebar>0.0:
            center=[0.2,info_height]
            if scalebar_pos is not None:
                center = scalebar_pos
            add_scalebar(ax,scalebar,center=center)
            ax.set_yticks([])
            ax.set_xticks([])
        return ims

    def plot(self,cbar_height=0.82,cbar_border=False,**kwargs):
        fig, ax = plt.subplots(1,1,figsize=(6,6),dpi=150)
        ims = self.plot_ax(ax,**kwargs)
        for i,k in enumerate(ims):
            # colorbar
            w = 0.65/len(ims)
            cax = fig.add_axes([0.1625+(0.05+w)*i, cbar_height, w, 0.02])
            cb = fig.colorbar(ims[k]["im"], cax=cax, orientation='horizontal')
            clabel = ims[k].get("field_display_name","")
            if ims[k]["units"]!="":
                clabel += " (%s)"%(r"$"+ims[k]["units"].replace("**","^")+"$")
            cblabelkwargs = {}
            if cbar_border:
                cb.outline.set_edgecolor('black')
                cb.ax.tick_params(color="white")
                cblabelkwargs = dict(color="white",fontsize=20)
            txt = cax.set_title(clabel,**cblabelkwargs)

            #cb.minorticks_off()
            if cbar_border:
                pass
                cb.ax.tick_params(labelsize=12,colors="white") 
                for l in cb.ax.get_xticklabels():
                    l.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground='black')])
                print(cb.ax.get_title())
                xlabel = cb.ax.set_title(cb.ax.get_title(),**cblabelkwargs)
                xlabel.set_path_effects([PathEffects.withStroke(linewidth=2.5, foreground='black')])


        return fig,ax

    def get_imarr(self,field):
        pnames = ["field"]+["field_%i"%i for i in range(1,6)]
        for i,pname in enumerate(pnames):
            try:
                if type(field)==tuple:
                    if self.params[pname]==field:
                        return i,self.imarrs["%i"%i]
                else:
                    if self.params[pname][1]==field:
                        return i,self.imarrs["%i"%i]
            except KeyError:
                pass
        print("could not find field '%s'"%field)
        return None,None

    def add_flamcontour(self,ax):
        extent = self.params["extent"]
        mid,imarr = self.get_imarr("flam")
        x = np.linspace(extent[0],extent[1],imarr.shape[0])
        y = np.linspace(extent[1],extent[2],imarr.shape[1])
        X, Y = np.meshgrid(x, y)
        CS = ax.contour(X, Y, imarr[::-1,:], levels = [0.1],colors='k')



    def save(self,fp):
        """Save PlotObj"""
        with h5py.File(fp, "w") as hf:
            grp = hf.create_group("imarrs")
            for i,k in enumerate(self.imarrs):
                grp.create_dataset(k, data=self.imarrs[k])
            for k in self.params:
                if self.params[k] is None:
                    continue
                hf.attrs[k] = self.params[k]
    def load(self,fp):
        """Load PlotObj"""
        with h5py.File(fp, "r") as hf:
            self.imarrs = {}
            grp = hf["imarrs"]
            for i,k in enumerate(grp):
                self.imarrs[k] = np.array(grp[k])
            for k in hf.attrs:
                self.params[k] = hf.attrs[k]
                if type(self.params[k])==np.ndarray:
                    # make those tuples (again), if they should have been arrays, 
                    # we had saved them as dataset not attribute!
                    self.params[k] = tuple(self.params[k])

def plot_time(ax,time,center=[0.85,0.15],color="white",textseparation=0.03,lw=4):
    txt = ax.text(center[0], center[1]-textseparation, "t=%.2f s"%time, fontsize=22,color=color,
                  fontname="DejaVu Sans",fontweight="bold", horizontalalignment='center',
                  verticalalignment='top',transform=ax.transAxes,zorder=10)
            #bbox=dict(boxstyle='square,pad=0.2', fc='black', ec='none',alpha=0.7))
    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='black')])

def slice_create(ds,field,center="c",width=(5000,"km"),debug=False,axis=None):
    """ Creates a slice and returns 2D numpy array of such."""
    resolution = 1600
    if axis is None: # TODO: Adjust for 1D/2D
        axis = yt.funcs.fix_axis("x",ds) # default for 3D
    if debug:
        resolution = 100 # lower resolution for faster debugging
    (bounds, center, display_center) = get_window_parameters(axis, center, width, ds)
    slc = ds.slice(axis, center[axis],center=center)
    frb = slc.to_frb(width, resolution)
    return frb[field]

def slice_plot(ds,field,width=(5000,"km"),zoffset=False,plotprops=None,debug=False):
    slc = yt.SlicePlot(ds,"x",field,width=width)
    if debug: # faster plotting; lower resolution as one measure for speed-up
        slc.set_buff_size((100,100)) #double resolution.
    else:
        slc.set_buff_size((1600,1600)) #double resolution.
    if zoffset:
        slc.set_center((0,width[0]/2),"km")
    
    # get plot properties for field if existent
    if plotprops is not None:
        fpp = dict(plotprops)
    else:
        fpp = fpps.get(field,{})
    
    # colormap 
    cmap = plt.matplotlib.cm.get_cmap(fpp.get("cmap","viridis"))
    so = fpp.get("cmap_set_over",None)
    if so:
        cmap.set_over(so)
    su = fpp.get("cmap_set_under",None)
    if su:
        cmap.set_under(su)
    slc.set_cmap(field,cmap)
    zlim = fpp.get("zlim",None)
    if zlim:
        slc.set_zlim(field,*zlim)
    zlinthresh = fpp.get("zlinthresh",None)
    slc.set_log(field, fpp.get("log",True),zlinthresh)
    return slc

def quad_plot(ds,fields,width=(5000,"km")):
    fig = plt.figure()
    
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (2, 2),
                    axes_pad = 0.9,
                    #label_mode = "1",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")
    
    #fields = [("flash",'soundspeed'), ("flash",'velocityz'), ("flash",'temperature'), ("flash",'density')]
    
    ps = []
    for field in fields:
        ps.append(slice_plot(ds,field,width=width))
    
    for i, field in enumerate(fields):
        plot = ps[i].plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        ps[i]._setup_plots()
    return fig

def triple_plot(ds,fields,lfpps={},**kwargs):
    fig = plt.figure(figsize=(15,5))
    
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (1, 3),
                    axes_pad = 0.9,
                    #label_mode = "1",
                    share_all = True,
                    cbar_location="right",
                    cbar_mode="each",
                    cbar_size="3%",
                    cbar_pad="0%")
    
    #fields = [("flash",'soundspeed'), ("flash",'velocityz'), ("flash",'temperature'), ("flash",'density')]
    
    ps = []
    for field in fields:
        fpp = lfpps.get(field,{})
        if not(fpp):
            fpp = fpps.get(field,{})
        ps.append(slice_plot(ds,field,plotprops=fpp,**kwargs))
    
    for i, field in enumerate(fields):
        plot = ps[i].plots[field]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
        ps[i]._setup_plots()
    return fig

def add_scalebar(ax,size,center=[0.20,0.15],color="white",textseparation=0.03,lw=4):
    axis_to_data = ax.transAxes + ax.transData.inverted()
    center_data = axis_to_data.transform(center)
    p = ax.plot([center_data[0]-size/2,center_data[0]+size/2],[center_data[1],center_data[1]],color=color,lw=lw,zorder=11)
    p[0].set_path_effects([PathEffects.withStroke(linewidth=6*lw/4.0, foreground='black')])
    txt = ax.text(center[0], center[1]-textseparation, "%i"%size+" km", fontsize=22,color=color,
                  fontname="DejaVu Sans",fontweight="bold", horizontalalignment='center',
                  verticalalignment='top',transform=ax.transAxes,zorder=10)
            #bbox=dict(boxstyle='square,pad=0.2', fc='black', ec='none',alpha=0.7))
    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='black')])


def get_norm(fpp):
    vmin,vmax = None,None
    if fpp.get("zlim",None) is not None:
        vmin = fpp["zlim"][0]
        vmax = fpp["zlim"][1]
    norm = Normalize(vmin=vmin,vmax=vmax)
    if fpp.get("log",False):
        norm = LogNorm(vmin=vmin,vmax=vmax)
    return norm

def get_cmap(fpp):
    cmap = plt.matplotlib.cm.get_cmap(fpp.get("cmap","viridis"))
    so = fpp.get("cmap_set_over",None)
    if so:
        cmap.set_over(so)
    su = fpp.get("cmap_set_under",None)
    if su:
        cmap.set_under(su)
    return cmap

