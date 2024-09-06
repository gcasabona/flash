import pandas as pd
import numpy as np

class flameprops(object):
    """Reads flame property files provided by Alexei Poludnenko"""
    def __init__(self,filepath):
        self.filepath = filepath
        self.fieldnames = []
        with open(filepath) as f:
            for l in f.readlines():
                if l.startswith("# ["):
                    idx = int(l.split("[")[1].split("]")[0]) # not really needed as ascending order anyway
                    fieldname = l.split("=")[1].strip()
                    self.fieldnames.append(fieldname)
        self.data = pd.read_csv(filepath,comment="#",names=self.fieldnames, delimiter=r"\s+",engine='python')
        self.data = self.data.sort_values(by=['density'])
    def get_value(self,density,field):
        return np.interp(density,self.data["density"],self.data[field])
    def get_flamespeed(self,density):
        return self.get_value(density,"laminar flame speed")
    
    
class flameprops_flash(object):
    """Reads flame property file in FLASH"""
    def __init__(self,filepath):
        self.filepath = filepath
        self.fieldnames = []
        with open(filepath) as f:
            l = f.readlines()[3].split()
            self.fieldnames = [s.replace('"',"") for s in l]
        self.data = pd.read_csv(filepath,skiprows=4,sep=" ",names=self.fieldnames, delimiter=r"\s+")
        #with open(filepath) as f:
        #    for l in f.readlines():
        #        print(l)
        #        if l.startswith("# ["):
         #           idx = int(l.split("[")[1].split("]")[0]) # not really needed as ascending order anyway
         #           fieldname = l.split("=")[1].strip()
        #            self.fieldnames.append(fieldname)
        #self.data = pd.read_csv(filepath,comment="#",names=self.fieldnames,delimiter="  ",engine='python')
        #self.data = self.data.sort_values(by=['density'])
    def get_value(self,density,field,x22=0.0,x12=0.5,logintp=True):
        mask = (self.data["X12"]==x12) & (self.data["X22"]==x22)
        if logintp:
            return 10**np.interp(np.log10(density),self.data[mask]['LOG(DENS)'],self.data[mask][field])
        else:
            return np.interp(density,10**self.data[mask]['LOG(DENS)'],10**self.data[mask][field])
    def get_flamespeed(self,density,x22=0.0,x12=0.5,logintp=True):
        return self.get_value(density,"LOG(S)",x22=x22,x12=x12,logintp=logintp)
