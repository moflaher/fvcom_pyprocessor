from __future__ import division, print_function
import numpy as np
import pyproj as pyp
import os

from load_mesh import load_functions
from save_fvcom import save_functions
from meshClass import MESH


class FVCOM:
    """
    An FVCOM object. Holds all information about mesh and input files.
    """
    
    def __init__(self, debug=False):
        """ Initialize FVCOM class """

        if debug: print('-Debug mode on-')
        self._debug = debug
        
        self._history = ['FVCOM object initialized.']
        self._savepath = ''
        
        
        # Define mesh loading functions
        self.load = load_functions(self)
        self.save = save_functions(self)
                
        
        return
    
    
    def set_savepath(self, savepath):
        """
        Specifies the savepath.
        """
        
        if self._debug: print('-Setting savepath-')
        
        if savepath.endswith('/'):
            self._savepath = savepath
        else:
            self._savepath = savepath + '/'  
                
        if not os.path.exists(self._savepath): os.makedirs(self._savepath)
        
        self._history.append('Savepath set')
        
        return
        
        
    def set_projection(self, projstr):
        """
        Specifies the projection and converts fields if available.
        """
        
        if not hasattr(self, 'mesh'):
            self.mesh = MESH
        
        if self._debug: print('-Setting projection.-')
        self.mesh._projstr = projstr
        try:
            self.mesh.proj = pyp.Proj(proj=projstr)
        except RuntimeError:
            print('Error in setting projection.')
        
        # Check for long/lat and x/y
        if hasattr(self.mesh, 'lon') and hasattr(self.mesh, 'lat'):
            ll = True
        else:
            ll = False
            
        if hasattr(self.mesh, 'x') and hasattr(self.mesh, 'y'):
            xy = True
        else:
            xy = False
        
        self._history.append('Projection set')
        
        # Done if have ll and xy
        if ll and xy:
            return
        
        # If ll and not xy set xy
        if ll and not xy:
            self.mesh.x, self.mesh.y = self.mesh.proj(self.mesh.lon, self.mesh.lat)
            self.mesh.nodexy = np.vstack([self.mesh.x, self.mesh.y]).T
            self._history.append('Projection used to convert coordinates')    
            
        # If xy and not ll set ll
        if not ll and xy:
            self.mesh.lon, self.mesh.lat = self.mesh.proj(self.mesh.x, self.mesh.y, inverse=True) 
            self.mesh.nodell = np.vstack([self.mesh.lon, self.mesh.lat]).T            
            self._history.append('Projection used to convert coordinates')        
        
        return
        

