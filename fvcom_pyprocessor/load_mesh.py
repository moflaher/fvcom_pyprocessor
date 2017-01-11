from __future__ import division, print_function
import numpy as np
import pyproj as pyp

from meshClass import MESH

class load_functions(object):
    """
    load class to define all of the loading functions.
    """
    
    def __init__(self,fv):      
        self.fv = fv  
        return

    def neifile(self, filename, isll=True):
        """
        Loads a .nei file. Assumes neifile is in long/lat.
        """  
        
        self = self.fv
        self.mesh = MESH
        
        try:
            if self._debug: print('-Opening nei file-')
            fp = open(filename, 'r')
        except IOError:
            print('Can not find ' + filename)
            return

        if self._debug: print('-Loading header-')
        self.mesh.nnodes = int(fp.readline())
        self.mesh.max_neighbours = int(fp.readline())
        self.mesh.bounding_box = np.array([float(x) for x in fp.readline().split()])
        
        if self._debug: print('-Loading data-')
        t_data = np.loadtxt(filename, skiprows=3, dtype='float64')
        fp.close()


        


        self.mesh._nodenumber = t_data[:,0].astype(int)
        self.mesh._boundary_code = t_data[:,3].astype(int)
        self.mesh._neighbours = t_data[:,5:].astype(int)
        
        self.mesh.h = t_data[:,4]

        if isll:
            if self._debug: print('-Using long/lat data-')
            self._history.append('Loaded lon/lat data')
            self.mesh.nodell = t_data[:,1:3]
            self.mesh.lon = t_data[:,1]
            self.mesh.lat = t_data[:,2]
        else:
            if self._debug: print('-Using x/y data-')
            self._history.append('Loaded x/y data')
            self.mesh.nodexy = t_data[:,1:3]
            self.mesh.x = t_data[:,1]
            self.mesh.y = t_data[:,2]
        
        if hasattr(self.mesh, 'proj'):
            if self._debug: print('-Using specified projection' +\
                                  'to convert coordinates-')
            self._history.append('Using specified projection to convert coordiantes')
            
            if isll:
                self.mesh.x, self.mesh.y = self.mesh.proj(self.mesh.lon, self.mesh.lat)
                self.mesh.nodexy = np.vstack([self.mesh.x, self.mesh.y]).T
            else:
                self.mesh.lon, self.mesh.lat = self.mesh.proj(self.mesh.x, self.mesh.y, inverse=True)
                self.mesh.nodell = np.vstack([self.mesh.lon, self.mesh.lat]).T
        else:        
            if isll:
                if self._debug: print('-No specified projection creating LCC ' +\
                                      'projection to convert coordinates to xy-')
                self._history.append('Created LCC projection to convert coordiantes')
                
                xmax=np.nanmax(self.mesh.lon)
                xmin=np.nanmin(self.mesh.lon)
                ymax=np.nanmax(self.mesh.lat)
                ymin=np.nanmin(self.mesh.lat)
                xavg = ( xmax + xmin ) * 0.5;
                yavg = ( ymax + ymin ) * 0.5;
                ylower = ( ymax - ymin ) * 0.25 + ymin;
                yupper = ( ymax - ymin ) * 0.75 + ymin;
                
                projstr='lcc +lon_0='+str(xavg)+' +lat_0='+str(yavg)+' +lat_1='+str(ylower)+' +lat_2='+str(yupper)
            
                self.set_projection(projstr)    
        
        self._history.append('Loaded neifile')
        
        return 
