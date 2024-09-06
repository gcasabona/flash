import yt
import numpy as np
from yt.units import dimensions as dims


# defining fields for dataset "ds" only as not always having these new data fields. Thus interactively import!


def _uptimescale(field,data):
    """time scale of op2mag change"""
    u = data.ds.units
    return data["op2mag"]/np.abs(data[("flash","updt")]*u.cm/u.s**2)
ds.add_field(('flash','uptimescale'), function=_uptimescale,units="s",
              display_name=r"$\frac{1}{u'}\left|\frac{du'}{dt}\right|$", sampling_type="cell",force_override=True)

def _wrinklingfactorTFI(field,data):
    """TFI wrinkling factor"""
    return data[("flash","fspd")]/data[("flash","slam")]
ds.add_field(('flash','wrinklingfactorTFI'), function=_wrinklingfactorTFI,
              display_name=r"$\Theta_\mathrm{TFI}$", sampling_type="cell",force_override=True)




