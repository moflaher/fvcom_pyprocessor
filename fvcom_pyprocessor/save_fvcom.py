from __future__ import division, print_function
import os
import numpy as np
import netCDF4 as n4
import time as ttime
from collections import OrderedDict as OD
import matplotlib.dates as dates

from utilities import get_nv


def _file_factory(self, casename, endname, filetype='ascii'):        
    """
    Code to create blank files.
    Takes a casename and endname. filetype is optional.
    """
    if casename==[]:
        casename = self.casename  
    filename = self._savepath + casename + endname            
    
    try:
        if self._debug: print('-Creating ' + filename + ' file-')
        if os.path.isfile(filename):
            print('File ' + filename + ' already exists.\nPlease remove file and rerun.')
            return None, filename
        if 'netcdf' in filetype:
            fp = n4.Dataset(filename, 'w',format='NETCDF3_CLASSIC')
        else:
            fp = open(filename, 'w')
        
        return fp, filename
    except IOError:
        print('Error creating ' + filename)
        return None, filename


class save_functions(object):
    """
    save class to define all of the saveing functions.
    """
    
    def __init__(self,fv):      
        self.fv = fv  
        return      

    
    def fvcom_grd(self, casename=[]):
        """
        Save an FVCOM_grd.dat file from an FVCOM object.
        """
        
        self = self.fv      
     
        if not hasattr(self.mesh, 'nv'):
            if self._debug: print('-Calculating nv-')
            self._history.append('Calculated nv')
            self.mesh.nv = get_nv(self.mesh)
            self.mesh.nele = len(self.mesh.nv)   
            
        fp, filename = _file_factory(self, casename, '_grd.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('Node Number = {}\n'.format(self.mesh.nnodes))
        fp.write('Cell Number = {}\n'.format(self.mesh.nele))
                
        for i in range(0,self.mesh.nele):
            fp.write('{} {} {} {} {}\n'.format(i+1,self.mesh.nv[i,0],
                                                   self.mesh.nv[i,1],
                                                   self.mesh.nv[i,2],0))

        for i in range(0,self.mesh.nnodes):
            fp.write('{} {} {} {}\n'.format(i+1,self.mesh.nodell[i,0],
                                                self.mesh.nodell[i,1],
                                                self.mesh.h[i]))
        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
    def fvcom_dep(self, casename=[], isll=False):
        """
        Save an FVCOM_dep.dat file from an FVCOM object.
        """
        self = self.fv
        
        fp, filename = _file_factory(self, casename, '_dep.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('Node Number = {}\n'.format(self.mesh.nnodes))

        if isll:
            for i in range(0,self.mesh.nnodes):
                fp.write('{} {} {} {}\n'.format(i+1,self.mesh.nodell[i,0],
                                                    self.mesh.nodell[i,1],
                                                    self.mesh.h[i]))
        else:
            for i in range(0,self.mesh.nnodes):
                fp.write('{} {} {} {}\n'.format(i+1,self.mesh.nodexy[i,0],
                                                    self.mesh.nodexy[i,1],
                                                    self.mesh.h[i]))
        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def fvcom_cor(self, casename=[], isll=False):
        """
        Save an FVCOM_cor.dat file from an FVCOM object.
        """
        self = self.fv
        
        fp, filename = _file_factory(self, casename, '_cor.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('Node Number = {}\n'.format(self.mesh.nnodes))

        if isll:
            for i in range(0,self.mesh.nnodes):
                fp.write('{} {} {}\n'.format(self.mesh.nodell[i,0],
                                             self.mesh.nodell[i,1],
                                             self.mesh.lat[i]))
        else:
            for i in range(0,self.mesh.nnodes):
                fp.write('{} {} {}\n'.format(self.mesh.nodexy[i,0],
                                             self.mesh.nodexy[i,1],
                                             self.mesh.lat[i]))
        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return


    def fvcom_spg(self, casename=[], spg=[], spgdist=[], spgval=[]):
        """
        Save an FVCOM_spg.dat file from an FVCOM object.
        Uses the specified spg values, or defaults to mesh.spg, then mesh.obc, then just writes empty spg file with a warning.
        """
        self = self.fv
        
        if spg==[]:
            if hasattr(self, 'spg'):
                spg = self.mesh.spg
            elif hasattr(self, 'obc'):
                spg = self.mesh.obc                
        if len(spg)==0:
            print('Warning saving empty _cage.dat file')
                
        # Error states
        if (len(spgdist) != len(spg)) or (len(spgval) != len(spg)):
            print('The length of spgdist, spgval, and spg must match.')
            return          
        
        fp, filename = _file_factory(self, casename, '_spg.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('Sponge Node Number = {}\n'.format(len(spg)))


        for i in range(0,len(spg)):
            fp.write('{} {} {}\n'.format(spg[i], spgdist[i], spgval[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
    def fvcom_obc(self, casename=[], obc=[]):
        """
        Save an FVCOM_obc.dat file from an FVCOM object.
        Uses the specified obc values, or defaults to mesh.obc then just writes empty obc file with a warning.
        """
        self = self.fv       
          
        if obc==[]:
            if hasattr(self, 'obc'):
                obc = self.mesh.obc           
        if len(obc)==0:
            print('Warning saving empty _obc.dat file')
                       
        
        fp, filename = _file_factory(self, casename, '_obc.dat')
        if fp == None: return
        
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        
        fp.write('OBC Node Number = {}\n'.format(len(obc)))

        for i in range(0,len(obc)):
            fp.write('{} {} {}\n'.format(i+1,obc[i],1))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return


    def fvcom_bfric(self, casename=[], bfric=[]):
        """
        Save an FVCOM_bfric.dat file from an FVCOM object.
        Uses the specified bfric values at each element.
        """
        self = self.fv   
         
        if bfric==[]:
            if hasattr(self, 'bfric'):
                bfric = self.mesh.bfric       
        
        if not ((len(bfric) == self.mesh.nnodes) or (len(bfric) == 0)):
            print('bfric must be the length of mesh.nnodes')
            return
                   
        fp, filename = _file_factory(self, casename, '_bfric.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        
        fp.write('BFRIC Node Number = {}\n'.format(self.mesh.nnodes))

        for i in range(0, len(bfric)):
            fp.write('{} {}\n'.format(i+1,bfric[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
        
    def fvcom_cage(self, casename=[], cage=[], cageval=[], cagedist=[]):
        """
        Save an FVCOM_cage.dat file from an FVCOM object.
        Uses the specified cage values, if empty writes empty cage file with a warning.
        Cage is define at the elements.
        Cagesdist is from the top down.
        """
        self = self.fv        
        
        if cage==[]:
            if hasattr(self, 'cage'):
                cage = self.mesh.cage           
        if len(cage)==0:
            print('Warning saving empty _cage.dat file')
                
        # Error states
        if (len(cagedist) != len(cage)) or (len(cageval) != len(cage)):
            print('The length of cagedist, cageval, and cage must match.')
            return            
            
        fp, filename = _file_factory(self, casename, '_cage.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('CAGE Node Number = {}\n'.format(len(cage)))


        for i in range(0,len(cage)):
            fp.write('{} {} {}\n'.format(cage[i], cageval[i], cagedist[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def fvcom_turbine(self, casename=[], turbine=[], turbineval=[], turbinedist=[]):
        """
        Save an FVCOM_turbine.dat file from an FVCOM object.
        Uses the specified turbine values, if empty writes empty turbine file with a warning.
        Turbine is defined at the elements.
        Turbinedist is from the bottom up.
        """
        self = self.fv        
        
        if turbine==[]:
            if hasattr(self, 'turbine'):
                turbine = self.mesh.turbine           
        if len(turbine)==0:
            print('Warning saving empty _turbine.dat file')
                
        # Error states
        if (len(turbinedist) is not len(turbine)) or (len(turbineval) is not len(turbine)):
            print('The length of turbinedist, turbineval, and turbine must match.')
            return
            
        fp, filename = _file_factory(self, casename, '_turbine.dat')
        if fp == None: return

            
        if self._debug: print('-Writing ' + filename + ' contents-')
        fp.write('TURBINE Node Number = {}\n'.format(len(turbine)))


        for i in range(0,len(turbine)):
            fp.write('{} {} {}\n'.format(turbine[i], turbineval[i], turbinedist[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def mesh_lon(self, casename=[], lon=[]):
        """
        Save a long.dat file from an FVCOM object.
        """
        self = self.fv        
        
        if lon==[]:
            if hasattr(self.mesh, 'lon'):
                lon = self.mesh.lon  
                
        if len(lon)==0:
            print('The length of lon must be non-zero.') 
            return
            
        if len(lon) != self.mesh.nnodes:
            print('Warning the length of lon does not equal self.mesh.nnodes') 

        fp, filename = _file_factory(self, casename, '_lon.dat')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        
        for i in range(0,len(lon)):
            fp.write('{}\n'.format(lon[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def mesh_lat(self, casename=[], lat=[]):
        """
        Save a lat.dat file from an FVCOM object.
        """
        self = self.fv  
        
        if lat==[]:
            if hasattr(self.mesh, 'lat'):
                lat = self.mesh.lat  
                
        if len(lat)==0:
            print('The length of lat must be non-zero') 
            return
            
        if len(lat) != self.mesh.nnodes:
            print('Warning the length of lat does not equal self.mesh.nnodes') 


        fp, filename = _file_factory(self, casename, '_lat.dat')
        if fp == None: return
        
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        
        for i in range(0,len(lat)):
            fp.write('{}\n'.format(lat[i]))

        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def neifile(self, casename=[]):
        """
        Save an nei-file file from an FVCOM object.
        """
        self = self.fv
        
        fp, filename = _file_factory(self, casename, '.nei')
        if fp == None: return
            
        if self._debug: print('-Writing ' + filename + ' contents-')
        
        fp.write('{:d}\n'.format(self.mesh.nnodes))
        fp.write('{:d}\n'.format(self.mesh.max_neighbours))
        fp.write('{:f} {:f} {:f} {:f}\n'.format(self.mesh.bounding_box[0],
                                                self.mesh.bounding_box[1], 
                                                self.mesh.bounding_box[2], 
                                                self.mesh.bounding_box[3]))
   

        for i in range(0,self.mesh.nnodes):
            fp.write('{:d} {:f} {:f} {:d} {:f} {}\n'.format(i+1, self.mesh.lon[i],
                                                            self.mesh.lat[i],
                                                            self.mesh._boundary_code[i],
                                                            self.mesh.h[i],
                                                            np.array_str(self.mesh._neighbours[i,])[1:-1]))
    
        fp.close()
        
        self._history.append('Wrote ' + filename)
        
        return
        
        
    def fvcom_spectide(self, casename=[]):
        """
        Save a FVCOM_spectide file from an FVCOM object.
        """
        self = self.fv
        
        if hasattr(self.mesh, 'spectide'):
            st = self.mesh.spectide
        else:
            print('mesh does not contain spectide dictionary')
            return
            
        if 'obc' not in st and hasattr(self.mesh, 'obc'):
            print('No spectide obc. Using mesh obc.')
            st['obc'] = self.mesh.obc
            
        # Make sure all fields that are required are in spectide
        req = ['components','obc','time_origin','tide_period',
               'tide_Ephase','tide_Eamp','equilibrium_tide_Eamp',
               'equilibrium_beta_love']
        checklist=[var for var in req if var not in st.keys() ]
        
        if len(checklist)>0:
            print('Cannot create spectide because the following required fields are missing:')
            print(checklist)
            return
        
        ncid, filename = _file_factory(self, casename, '_spectide.nc', filetype='netcdf')
        if ncid == None: return
            
        ncattrs = OD()    
        ncattrs['type'] = 'FVCOM SPECTRAL ELEVATION FORCING FILE'
        ncattrs['title'] = 'Spectral tidal boundary input'
        ncattrs['components'] = ''
        ncattrs['history'] = 'File created using fvcom_spectide from the Python fvcom_pyprocessor toolbox on ' + ttime.ctime()
                         
        components = ''
        for component in st['components']:
            components += component + ' '
        ncattrs['components'] = components.strip()
        if 'title' in st:
            ncattrs['title'] = st['title']
            
        # Set global ncfile attributes    
        for attname in ncattrs:
            setattr(ncid,attname,ncattrs[attname])
            
        ncdim = OD()
        ncdim['nobc'] = len(st['obc'])
        ncdim['tidal_components'] = len(st['components'])
        ncdim['TideLen'] = 26
        ncdim['DateStrLen'] = 19

        
        # Set the ncfile dimensions
        for dimname in ncdim:
            ncid.createDimension(dimname,ncdim[dimname])
        
        
        # Create variables and set their attributes
        ncvars = {}
        ncvars['obc'] = ncid.createVariable('obc_nodes', 'i4', ('nobc'))
        ncvars['obc'].setncattr('long_name', 'Open Boundary Node Number')
        ncvars['obc'].setncattr('grid', 'obc_grid')

        ncvars['tide_period'] = ncid.createVariable('tide_period', 'f4', ('tidal_components'))
        ncvars['tide_period'].setncattr('long_name', 'Tide Angular Period')
        ncvars['tide_period'].setncattr('units', 'seconds')
        
        ncvars['tide_Eref'] = ncid.createVariable('tide_Eref', 'f4', ('nobc'))
        ncvars['tide_Eref'].setncattr('long_name', 'Tidal Elevation Reference Level')
        ncvars['tide_Eref'].setncattr('units', 'meters')
        
        ncvars['tide_Ephase'] = ncid.createVariable('tide_Ephase', 'f4', ('nobc','tidal_components'))
        ncvars['tide_Ephase'].setncattr('long_name', 'Tidal Elevation Phase Angle')
        ncvars['tide_Ephase'].setncattr('units', 'degrees, time of maximum elevation with respect to chosen time origin')
        
        ncvars['tide_Eamp'] = ncid.createVariable('tide_Eamp', 'f4', ('nobc','tidal_components'))
        ncvars['tide_Eamp'].setncattr('long_name', 'Tidal Elevation Amplitude')
        ncvars['tide_Eamp'].setncattr('units', 'meters')
        
        ncvars['equilibrium_tide_Eamp'] = ncid.createVariable('equilibrium_tide_Eamp', 'f4', ('tidal_components'))
        ncvars['equilibrium_tide_Eamp'].setncattr('long_name', 'Equilibrium Tidal Elevation Amplitude')
        ncvars['equilibrium_tide_Eamp'].setncattr('units', 'meters')
        
        ncvars['equilibrium_beta_love'] = ncid.createVariable('equilibrium_beta_love', 'f4', ('tidal_components'))
        ncvars['equilibrium_beta_love'].setncattr('formula', 'beta=1+klove-hlvoe')
        
        ncvars['equilibrium_tide_type'] = ncid.createVariable('equilibrium_tide_type', 'c', ('tidal_components','TideLen'))
        ncvars['equilibrium_tide_type'].setncattr('long_name', 'formula')
        ncvars['equilibrium_tide_type'].setncattr('units', 'beta=1+klove-hlove')
        
        # Create time variable
        ncvars['time_origin'] = ncid.createVariable('time_origin', 'c', ('DateStrLen'))
        ncvars['time_origin'].setncattr('long_name', 'time')
        ncvars['time_origin'].setncattr('units', 'yyyy-mm-dd HH:MM:SS')
        ncvars['time_origin'].setncattr('time_zone', 'UTC')    
        ncvars['time_origin'].setncattr('comments', 'Tidal Harmonic origin_date')
        
        
        # Fill variables
        st['equilibrium_tide_type'] = np.empty((len(st['tide_period']),26), dtype=str)

        for i,period in enumerate(st['tide_period']):
            if period <= 13*3600:
                st['equilibrium_tide_type'][i,] = list('SEMIDIURNAL'.ljust(26))
            elif (period > 13*3600) and (period < 28*3600):
                st['equilibrium_tide_type'][i,] = list('DIURNAL'.ljust(26))
            else:
                print('A period has been specified that is unsupported by FVCOM')
                st['equilibrium_tide_type'][i,] = list('LONG PERIOD'.ljust(26))
                
        for key in ncvars:
            ncvars[key][:] = st[key][:]
            
        ncid.close()     
                    
        return
        

    def fvcom_initTS(self, casename=[]):
        """
        Save a FVCOM initial TS file from an FVCOM object.
        """
        self = self.fv
        
        if hasattr(self.mesh, 'initTS'):
            ts = self.mesh.initTS
        else:
            print('mesh does not contain initTS dictionary')
            return
                       
        # Make sure all fields that are required are in spectide
        req = ['Times','zsl','tsl','ssl']
        checklist=[var for var in req if var not in ts.keys() ]
        
        if len(checklist)>0:
            print('Cannot create initTS because the following required fields are missing:')
            print(checklist)
            return
        
        ncid, filename = _file_factory(self, casename, '_initTS.nc', filetype='netcdf')
        if ncid == None: return
            
        ncattrs = OD()    
        ncattrs['type'] = 'FVCOM TEMP/SAL INITIALIZATION FILE'
        ncattrs['title'] = 'temp/sal initialiazation values'
        ncattrs['history'] = 'File created using fvcom_initTS from the Python fvcom_pyprocessor toolbox on ' + ttime.ctime()
            
        if 'title' in ts:
            ncattrs['title'] = ts['title']
            
        # Set global ncfile attributes    
        for attname in ncattrs:
            setattr(ncid,attname,ncattrs[attname])
            
        ncdim = OD()
        ncdim['DateStrLen'] = 19
        ncdim['time'] = None
        ncdim['node'] = self.mesh.nnodes
        ncdim['ksl'] = len(ts['zsl'])
        
        # Set the ncfile dimensions
        for dimname in ncdim:
            ncid.createDimension(dimname,ncdim[dimname])
        
        
        # Create variables and set their attributes
        ncvars = {}
        ncvars['Times'] = ncid.createVariable('Times', 'c', ('time','DateStrLen'))
        ncvars['Times'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['Times'].setncattr('time_zone', 'UTC')

        ncvars['time'] = ncid.createVariable('time', 'f4', ('time'))
        ncvars['time'].setncattr('long_name', 'Time')
        ncvars['time'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['time'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime'] = ncid.createVariable('Itime', 'i4', ('time'))
        ncvars['Itime'].setncattr('long_name', 'Time')
        ncvars['Itime'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['Itime'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime2'] = ncid.createVariable('Itime2', 'i4', ('time'))
        ncvars['Itime2'].setncattr('long_name', 'Time')
        ncvars['Itime2'].setncattr('units', 'msec since 00:00:00')
        ncvars['Itime2'].setncattr('time_zone', 'UTC')
        
        ncvars['zsl'] = ncid.createVariable('zsl', 'f4', ('ksl'))
        ncvars['zsl'].setncattr('long_name', 'Standard Z Levels, Positive Up')
        ncvars['zsl'].setncattr('units', 'm')

        ncvars['tsl'] = ncid.createVariable('tsl', 'f4', ('time','ksl','node'))
        ncvars['tsl'].setncattr('long_name', 'Observed Temperature Profiles')
        ncvars['tsl'].setncattr('units', 'degrees C')
        
        ncvars['ssl'] = ncid.createVariable('ssl', 'f4', ('time','ksl','node'))
        ncvars['ssl'].setncattr('long_name', 'Observed Salinity Profiles')
        ncvars['ssl'].setncattr('units', 'PSU')
        
        
        # Fill variables
        refdate = dates.datestr2num('1858-11-17 00:00:00')
        date = dates.datestr2num(ts['Times'])
        
        ts['time'] = date - refdate
        ts['Itime'] = np.floor(ts['time'])
        ts['Itime2'] = np.array([0])
        
        Times = np.empty((len(ts['Times']),19),dtype=str)        
        for i,time in enumerate(ts['Times']):
            Times[i,:] = list(time)
        ts['Times'] = Times        
        
        
        for key in ncvars:
            ncvars[key][:] = ts[key][:]
            
        ncid.close()     
                    
        return
        
        
    def fvcom_wind(self, casename=[]):
        """
        Save a FVCOM wind forcing file from an FVCOM object.
        """
        self = self.fv
        
        if hasattr(self.mesh, 'wind'):
            wnd = self.mesh.wind
        else:
            print('mesh does not contain wind dictionary')
            return
                       
        # Make sure all fields that are required are in spectide
        req = ['Times','U10','V10']
        checklist=[var for var in req if var not in wnd.keys() ]
        
        if len(checklist)>0:
            print('Cannot create wind because the following required fields are missing:')
            print(checklist)
            return
        
        ncid, filename = _file_factory(self, casename, '_wind.nc', filetype='netcdf')
        if ncid == None: return
            
        ncattrs = OD()    
        ncattrs['type'] = 'FVCOM WIND FORCING FILE'
        ncattrs['title'] = 'wind forcing values'
        ncattrs['history'] = 'File created using fvcom_wind from the Python fvcom_pyprocessor toolbox on ' + ttime.ctime()
            
        if 'title' in wnd:
            ncattrs['title'] = wnd['title']
            
        # Set global ncfile attributes    
        for attname in ncattrs:
            setattr(ncid,attname,ncattrs[attname])
            
        ncdim = OD()
        ncdim['time'] = None
        ncdim['nele'] = self.mesh.nele
        ncdim['DateStrLen'] = 19
        
        # Set the ncfile dimensions
        for dimname in ncdim:
            ncid.createDimension(dimname,ncdim[dimname])
        
        
        # Create variables and set their attributes
        ncvars = {}
        ncvars['Times'] = ncid.createVariable('Times', 'c', ('time','DateStrLen'))
        ncvars['Times'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['Times'].setncattr('time_zone', 'UTC')

        ncvars['time'] = ncid.createVariable('time', 'f4', ('time'))
        ncvars['time'].setncattr('long_name', 'Time')
        ncvars['time'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['time'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime'] = ncid.createVariable('Itime', 'i4', ('time'))
        ncvars['Itime'].setncattr('long_name', 'Time')
        ncvars['Itime'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['Itime'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime2'] = ncid.createVariable('Itime2', 'i4', ('time'))
        ncvars['Itime2'].setncattr('long_name', 'Time')
        ncvars['Itime2'].setncattr('units', 'msec since 00:00:00')
        ncvars['Itime2'].setncattr('time_zone', 'UTC')
        
        ncvars['U10'] = ncid.createVariable('U10', 'f4', ('time','nele'))
        ncvars['U10'].setncattr('long_name', 'Eastward Wind speed at 10m height')
        ncvars['U10'].setncattr('standard_name', 'Wind speed')
        ncvars['U10'].setncattr('units', 'm/s')
        ncvars['U10'].setncattr('type', 'data')

        ncvars['V10'] = ncid.createVariable('V10', 'f4', ('time','nele'))
        ncvars['V10'].setncattr('long_name', 'Northward Wind speed at 10m height')
        ncvars['V10'].setncattr('standard_name', 'Wind speed')
        ncvars['V10'].setncattr('units', 'm/s')
        ncvars['V10'].setncattr('type', 'data')
        
        
        # Fill variables
        refdate = dates.datestr2num('1858-11-17 00:00:00')
        date = dates.datestr2num(wnd['Times'])
        
        wnd['time'] = date - refdate
        wnd['Itime'] = np.floor(wnd['time'])
        wnd['Itime2'] = np.array([0])
        
        Times = np.empty((len(wnd['Times']),19),dtype=str)        
        for i,time in enumerate(wnd['Times']):
            Times[i,:] = list(time)
        wnd['Times'] = Times    
        
        for key in ncvars:
            ncvars[key][:] = wnd[key][:]
            
        ncid.close()     
                    
        return
        
        
    def fvcom_elevtide(self, casename=[]):
        """
        Save a FVCOM elevation boundary forcing file from an FVCOM object.
        """
        self = self.fv
        
        if hasattr(self.mesh, 'elevtide'):
            eld = self.mesh.elevtide
        else:
            print('mesh does not contain elevtide dictionary')
            return
                       
        # Make sure all fields that are required are in elevtide
        req = ['nobc','Times','elevation','obc_nodes']
        checklist=[var for var in req if var not in eld.keys() ]
        
        if len(checklist)>0:
            print('Cannot create elevation forcing because the following required fields are missing:')
            print(checklist)
            return
        
        ncid, filename = _file_factory(self, casename, '_elevtide.nc', filetype='netcdf')
        if ncid == None: return
            
        ncattrs = OD()    
        ncattrs['type'] = 'FVCOM TIME SERIES ELEVATION FORCING FILE'
        ncattrs['title'] = 'surface elevation boundary input'
        ncattrs['history'] = 'File created using fvcom_elevtide from the Python fvcom_pyprocessor toolbox on ' + ttime.ctime()
            
        if 'title' in eld:
            ncattrs['title'] = eld['title']
            
        # Set global ncfile attributes    
        for attname in ncattrs:
            setattr(ncid,attname,ncattrs[attname])
            
        ncdim = OD()
        ncdim['DateStrLen'] = 26
        ncdim['time'] = None
        ncdim['nobc'] = eld['nobc']
        
        # Set the ncfile dimensions
        for dimname in ncdim:
            ncid.createDimension(dimname,ncdim[dimname])
        
        
        # Create variables and set their attributes
        ncvars = {}
        ncvars['Times'] = ncid.createVariable('Times', 'c', ('time','DateStrLen'))
        ncvars['Times'].setncattr('time_zone', 'UTC')

        ncvars['time'] = ncid.createVariable('time', 'f4', ('time'))
        ncvars['time'].setncattr('long_name', 'Time')
        ncvars['time'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['time'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['time'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime'] = ncid.createVariable('Itime', 'i4', ('time'))
        ncvars['Itime'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['Itime'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['Itime'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime2'] = ncid.createVariable('Itime2', 'i4', ('time'))
        ncvars['Itime2'].setncattr('units', 'msec since 00:00:00')
        ncvars['Itime2'].setncattr('time_zone', 'UTC')
        
        ncvars['elevation'] = ncid.createVariable('elevation', 'f4', ('time','nobc'))
        ncvars['elevation'].setncattr('long_name', 'Open Boundary Elevation')
        ncvars['elevation'].setncattr('units', 'meters')
        
        ncvars['obc_nodes'] = ncid.createVariable('obc_nodes', 'i4', ('nobc'))
        ncvars['obc_nodes'].setncattr('long_name', 'Open Boundary Node Number')
        ncvars['obc_nodes'].setncattr('grid', 'obc_grid')

        ncvars['iint'] = ncid.createVariable('iint', 'i4', ('time'))
        ncvars['iint'].setncattr('long_name', 'internal mode iteration number')
        
        
        # Fill variables
        refdate = dates.datestr2num('1858-11-17 00:00:00')
        date = dates.datestr2num(eld['Times'])
        
        eld['time'] = date - refdate
        eld['Itime'] = np.floor(eld['time'])
        eld['Itime2'] = (eld['time'] - eld['Itime']) * 60 * 60 * 24 * 1000
        eld['iint'] = np.arange(1,len(eld['time'])+1,dtype=int)         
        
        Times = np.empty((len(eld['Times']),26),dtype=str)        
        for i,time in enumerate(eld['Times']):
            tlist = list(time)
            if len(tlist)==26:
                Times[i,:] = tlist
            elif len(tlist)>26:
                print('Times value has more than 26 characters: truncating')
                Times[i,:] = tlist[:26]
            elif len(tlist)<26:
                print('Times value has less than 26 characters: padding')
                Times[i,:] = tlist[:26] + '0' * (26-len(tlist))
            else:
                print('Something is very wrong with a Times value')
                
        eld['Times'] = Times       
        
                
        for key in ncvars:
            ncvars[key][:] = eld[key][:]
            
        ncid.close()     
                    
        return
        
        
    def fvcom_tsobc(self, casename=[]):
        """
        Save a FVCOM temp/salanity open boundary forcing file from an FVCOM object.
        """
        self = self.fv
        
        if hasattr(self.mesh, 'tsobc'):
            tsobc = self.mesh.tsobc
        else:
            print('mesh does not contain tsobc dictionary')
            return
                       
        # Make sure all fields that are required are in tsobc
        req = ['nobc','obc_nodes','Times','obc_temp','obc_salinity','siglay','siglev']
        checklist=[var for var in req if var not in tsobc.keys() ]
        
        if len(checklist)>0:
            print('Cannot create tsobc forcing because the following required fields are missing:')
            print(checklist)
            return
        
        ncid, filename = _file_factory(self, casename, '_tsobc.nc', filetype='netcdf')
        if ncid == None: return
            
        ncattrs = OD()    
        ncattrs['type'] = 'FVCOM TIME SERIES OBS TS FORCING FILE'
        ncattrs['title'] = 'FVCOM HYDROGRAPHIC OPEN BOUNDARY FORCING FILE'
        ncattrs['history'] = 'File created using fvcom_tsobc from the Python fvcom_pyprocessor toolbox on ' + ttime.ctime()
            
        if 'title' in tsobc:
            ncattrs['title'] = tsobc['title']
            
        # Set global ncfile attributes    
        for attname in ncattrs:
            setattr(ncid,attname,ncattrs[attname])
            
        ncdim = OD()
        ncdim['DateStrLen'] = 26
        ncdim['time'] = None
        ncdim['nobc'] = tsobc['nobc']
        ncdim['siglay'] = len(tsobc['siglay'])
        ncdim['siglev'] = len(tsobc['siglay']) + 1
        
        # Set the ncfile dimensions
        for dimname in ncdim:
            ncid.createDimension(dimname,ncdim[dimname])
        
        
        # Create variables and set their attributes
        ncvars = {}
        ncvars['Times'] = ncid.createVariable('Times', 'c', ('time','DateStrLen'))
        ncvars['Times'].setncattr('time_zone', 'UTC')

        ncvars['time'] = ncid.createVariable('time', 'f4', ('time'))
        ncvars['time'].setncattr('long_name', 'Time')
        ncvars['time'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['time'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['time'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime'] = ncid.createVariable('Itime', 'i4', ('time'))
        ncvars['Itime'].setncattr('units', 'days since 1858-11-17 00:00:00')
        ncvars['Itime'].setncattr('format', 'Modified Julian Day (MJD)')
        ncvars['Itime'].setncattr('time_zone', 'UTC')
        
        ncvars['Itime2'] = ncid.createVariable('Itime2', 'i4', ('time'))
        ncvars['Itime2'].setncattr('units', 'msec since 00:00:00')
        ncvars['Itime2'].setncattr('time_zone', 'UTC')
                
        ncvars['obc_nodes'] = ncid.createVariable('obc_nodes', 'i4', ('nobc'))
        ncvars['obc_nodes'].setncattr('long_name', 'Open Boundary Node Number')
        ncvars['obc_nodes'].setncattr('grid', 'obc_grid')
        
        ncvars['obc_h'] = ncid.createVariable('obc_h', 'f4', ('nobc'))
        ncvars['obc_h'].setncattr('long_name', 'Open Boundary Node Depth')
        ncvars['obc_h'].setncattr('units', 'm')
        ncvars['obc_h'].setncattr('grid', 'obc_grid')
        
        ncvars['obc_siglev'] = ncid.createVariable('obc_siglev', 'f4', ('siglev','nobc'))
        ncvars['obc_siglev'].setncattr('long_name', 'ocean_sigma/general_coordinate')
        ncvars['obc_siglev'].setncattr('grid', 'obc_grid')
        
        ncvars['obc_siglay'] = ncid.createVariable('obc_siglay', 'f4', ('siglay','nobc'))
        ncvars['obc_siglay'].setncattr('long_name', 'ocean_sigma/general_coordinate')
        ncvars['obc_siglay'].setncattr('grid', 'obc_grid')
        
        ncvars['obc_temp'] = ncid.createVariable('obc_temp', 'f4', ('time','siglay','nobc'))
        ncvars['obc_temp'].setncattr('long_name', 'Sea Water Temperature')
        ncvars['obc_temp'].setncattr('units', 'Celsius')
        ncvars['obc_temp'].setncattr('grid', 'obc_grid')

        ncvars['obc_salinity'] = ncid.createVariable('obc_salinity', 'f4', ('time','siglay','nobc'))
        ncvars['obc_salinity'].setncattr('long_name', 'Sea Water Salinity')
        ncvars['obc_salinity'].setncattr('units', 'PSU')
        ncvars['obc_salinity'].setncattr('grid', 'obc_grid')


        # Fill variables
        refdate = dates.datestr2num('1858-11-17 00:00:00')
        date = dates.datestr2num(tsobc['Times'])
        
        tsobc['time'] = date - refdate
        tsobc['Itime'] = np.floor(tsobc['time'])
        tsobc['Itime2'] = (tsobc['time'] - tsobc['Itime']) * 60 * 60 * 24 * 1000        
        
        Times = np.empty((len(tsobc['Times']),26),dtype=str)        
        for i,time in enumerate(tsobc['Times']):
            tlist = list(time)
            if len(tlist)==26:
                Times[i,:] = tlist
            elif len(tlist)>26:
                print('Times value has more than 26 characters: truncating')
                Times[i,:] = tlist[:26]
            elif len(tlist)<26:
                print('Times value has less than 26 characters: padding')
                Times[i,:] = tlist[:26] + '0' * (26-len(tlist))
            else:
                print('Something is very wrong with a Times value')
                
        tsobc['Times'] = Times 
        
        tsobc['obc_h'] = self.mesh.h[tsobc['obc_nodes'] - 1]            
        tsobc['obc_siglay'] = np.tile(tsobc['siglay'], (tsobc['nobc'],1)).T
        tsobc['obc_siglev'] = np.tile(tsobc['siglev'], (tsobc['nobc'],1)).T
                
        for key in ncvars:
            ncvars[key][:] = tsobc[key][:]
            
        ncid.close()     
                    
        return
