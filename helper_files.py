import configparser
import os

class runorganizer(object):
    """Minimal utility to fetch runs and relevant files to drag them into yt."""
    def __init__(self,configpath):
        self.config = configparser.ConfigParser()
        cfl = self.config.read(configpath)
        if len(cfl)==0:
            raise FileNotFoundError("Cannot find config.")
        self.update_runlist()
    def update_runlist(self):
        path = self.config.get("DEFAULT","runpath")
        folders = [os.path.join(path, o) for o in os.listdir(path) if os.path.isdir(os.path.join(path,o))]
        runs = [f for f in folders if os.path.isfile(os.path.join(f,"flash.par"))]
        self.runpaths = {r.split("/")[-1]:r for r in runs}
        return list(self.runpaths.keys())
