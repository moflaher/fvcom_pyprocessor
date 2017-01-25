from __future__ import division, print_function
import numpy as np
import xarray as xr
import numpy.ma as ma
from scipy.interpolate import LinearNDInterpolator as lndi
import matplotlib.dates as dates
import copy

def get_nv(self):
    """
    Finds the element vertex's.
    """
    
    
    neighbours_orig = self._neighbours
    nnodes = self.nnodes
    maxnei = self.max_neighbours
    

    try:
        import pyximport; pyximport.install()
        import get_nv as gnv
        
        nv = gnv.get_nvc(neighbours_orig, nnodes, maxnei)

    except:
        print('There was an issue with during using cython falling back to python.')
    
        nv = np.empty((len(neighbours_orig)*2, 3))    
        
        neighbours = copy.deepcopy(neighbours_orig)

        kk=0
        for i in range(nnodes-2):
            nei_cnt = 1
            for ii in range(maxnei-1):
                if neighbours[i, ii+1]==0:
                    break
                nei_cnt = ii+1    
                if neighbours[i, ii] <= (i+1):
                    continue
                if neighbours[i, ii+1] <= (i+1):
                    continue   
                for j in range(maxnei):
                    if neighbours[neighbours[i, ii]-1, j] != neighbours[i, ii+1]:
                        continue
                    nv[kk, :] = [i+1, neighbours[i,ii], neighbours[i, ii+1]]
                    kk = kk+1
                    break

            if (nei_cnt > 1):
                for j in range(maxnei):
                    if neighbours[i, 0] <= (i+1):
                        break
                    if neighbours[i, nei_cnt] <= (i+1):
                        break
                    if neighbours[neighbours[i, 0]-1, j] == 0:
                        break    
                    if neighbours[neighbours[i, 0]-1, j] != neighbours[i, nei_cnt]:
                        continue
                    nv[kk, :] = [i+1, neighbours[i, nei_cnt], neighbours[i, 0]]
                    kk = kk+1
                    break
                     
        nv = np.delete(nv, np.s_[kk:], axis=0)
        nv = (nv-1).astype(int)   
                
    return nv
    
    
def load_nemo(maskname, coordsname, outputname, varlist=[]):
    """
    Loads a nemo coordinate file and output file.
    """
    
    nemo = {}
    nemo['coords']={}
    nemo['mask']={}
    nemo['rawdata']={}
    nemo['mask_filename']=maskname
    nemo['coordinate_filename']=coordsname
    nemo['output_filename']=outputname
    
    # Load mask with xr netcdf
    ncid = xr.open_dataset(nemo['mask_filename'])
    for key in ncid.variables.keys():
        nemo['mask'][key]=ncid.variables[key]
    
    # Load coordinates with xr netcdf
    ncid = xr.open_dataset(nemo['coordinate_filename'])
    for key in ncid.variables.keys():
        nemo['coords'][key]=ncid.variables[key]
        
    # Load output with scipy netcdf
    ncid = xr.open_dataset(nemo['output_filename'])
    for key in ncid.variables.keys():
        nemo['rawdata'][key]=ncid.variables[key]
        
    # Set convenience fields
    nemo['LON']=nemo['coords']['nav_lon']
    nemo['LAT']=nemo['coords']['nav_lat']
    nemo['lon']=np.ravel(nemo['LON'])
    nemo['lat']=np.ravel(nemo['LAT'])
    nemo['x'],nemo['y'],nemo['proj']=lcc(nemo['lon'],nemo['lat'])
    nemo['xy']=np.array([nemo['x'],nemo['y']]).T
    nemo['ll']=np.array([nemo['lon'],nemo['lat']]).T
        
    # Don't set mask for all data. Causes large ram usage.
    # Not sure if possible to have memory mapped masks....
    
    
    return nemo
    
    
def interp_nemo(nemo, points, fieldname, time, layer=None, mymask=None):
    """
    Interpolation nemo data field to a set of points.
    """
    
    # If no specified mask try and pick mask based on fieldname.
    if mymask is None:
        masks = makemasks(nemo)
        mymask = masks[fieldname]
        
        
    # Sigh... slicing multi-dim arrays is not working the same as ipython...
    # Grrrrrrr....
    # So if we have layers then slice layer and time together else just time...
    
    # Then slice layer if specified
    if layer is not None:
        tdata = nemo['rawdata'][fieldname][time,layer,]
        mymask = ~mymask[layer,:].astype(bool)
    else:
        tdata = nemo['rawdata'][fieldname][time,]
            
    datain = ma.masked_array(tdata, mask=mymask, fill_value=np.nan)
    interpf = lndi(nemo['ll'], np.ravel(datain), fill_value=np.nan)
    values = interpf(points[:,0], points[:,1])
          
    return values
    
    
def lcc(lon,lat):
    """
    Given a lon lat converts to x,y and return them and the projection     
    """
    
    try:
        import pyproj as pyp
    except ImportError:
        print("pyproj is not installed, please install pyproj.")
        return
    
    
    # Define the lcc projection
    xmax = np.nanmax(lon)
    xmin = np.nanmin(lon)
    ymax = np.nanmax(lat)
    ymin = np.nanmin(lat)
    xavg = ( xmax + xmin ) * 0.5;
    yavg = ( ymax + ymin ) * 0.5;
    ylower = ( ymax - ymin ) * 0.25 + ymin;
    yupper = ( ymax - ymin ) * 0.75 + ymin;
    
    projstr = 'lcc +lon_0='+str(xavg)+' +lat_0='+str(yavg)+' +lat_1='+str(ylower)+' +lat_2='+str(yupper)
    proj = pyp.Proj(proj=projstr)
    
    x, y = proj(lon,lat)     
    
    return x, y, proj
    

def makemasks(nemo):
    """
    Create dictionary of masks.
    """
    
    masks={}
    
    try:
        masks['votemper']=nemo['mask']['tmask'][0,]
    except KeyError:
        print('Missing tmask, no mask set for votemper.')
    try:
        masks['vosaline']=nemo['mask']['tmask'][0,]
    except KeyError:
        print('Missing tmask, no mask set for vosaline.')
    try:
        # This is a hair hacky should find real mask
        a = np.min(nemo['rawdata']['sossheig'],axis=0)
        b = np.max(nemo['rawdata']['sossheig'],axis=0)
        masks['sossheig'] = a==b
    except KeyError:
        print('Can''t create sossheig mask, no mask set for sossheig.')
    
    
    return masks
    
    
def find_depth(nemo, depth):
    """
    Find the layer that is closest to each depth
    """
    
    h = nemo['rawdata']['deptht']
    idx = np.array([],dtype=np.int)
    for val in depth:
        tidx = np.argmin(np.fabs(h-val))
        idx=np.append(idx,tidx)
        
    return idx   
    
    
def find_time(nemo, time):
    """
    Find the time that is closest to each time
    """
    
    if type(time) is not np.ndarray:
        time = np.array(time, dtype=np.datetime64)
    if type(time.flat[0]) is not np.datetime64:
        time = time.astype(np.datetime64)

    times = nemo['rawdata']['time_counter']
    idx = np.array([],dtype=np.int)
    for val in time:
        tidx = np.argmin(np.abs(times-val))
        idx=np.append(idx,tidx)
        
    return idx  
    
def fix_nanzero(mesh, datain, fixnan=True, fixzero=True, cutoff=0):
    """
    Given a data set on a grid spread data to nans and zeros
    Note: only works with node data currently.
    """
    
    data = copy.deepcopy(datain)
    
    def idx_nanzero(data, fixnan, fixzero):
        idx=np.array([],dtype=int)
        if fixnan==True:
            idx = np.append(idx,np.argwhere(np.isnan(data)==1))
        if fixzero==True:
            idx = np.append(idx,np.argwhere(data<cutoff))
        return idx
    
    idx = idx_nanzero(data, fixnan, fixzero)
        
    while len(idx)>0:
        for i in idx:
            neis = (mesh._neighbours[i,]-1).astype(int)
            neis = neis[neis!=-1]            
            nodes = data[neis]
            nidx = idx_nanzero(nodes, fixnan, fixzero)
            bidx = np.ones(neis.shape,dtype=bool)
            bidx[nidx]=False
            if np.sum(bidx)>0:
                data[i] = np.mean(nodes[bidx]) 
        idx = idx_nanzero(data, fixnan, fixzero)
        
    return data  
    

def sort_boundary(mesh):
    """
    Given a mesh it sorts the external boundary and returns it as an array.
    """
    
    boundary = 1
    
    bcode = mesh._boundary_code
    nodenumber = mesh._nodenumber
    neighbours = mesh._neighbours
    
    nn = copy.deepcopy(nodenumber[bcode==boundary]).astype(int)
    nnei = copy.deepcopy(neighbours[bcode==boundary]).astype(int)
    
    #find the neighbour of the first node
    idx = np.argwhere(nnei==nn[0])[0][0]
    
    #have to use temp values with copy as the standard swap doesn't work when things are swapped again and again.
    #there must be a more python way to hand that....
    tmpval = nn[1].copy()
    nn[1] = nn[idx]
    nn[idx] = tmpval    
    tmpval = nnei[1,:].copy()
    nnei[1,:] = nnei[idx,:]
    nnei[idx,:] = tmpval
   
    for i in range(1,len(nn)-1):
        for j in range(mesh.max_neighbours):
            nei = nnei[i,j]
            if nei==0: continue
            idx = np.argwhere(nn[(i+1):]==nei)

            if len(idx)==1:
                tmpval = nn[(i+1)].copy()
                nn[(i+1)] = nn[(idx+i+1)]
                nn[(idx+i+1)] = tmpval                
                tmpval = nnei[(i+1),:].copy()
                nnei[(i+1),:] = nnei[(idx+i+1),:]
                nnei[(idx+i+1),:] = tmpval
                break
               
    return nn
    
    

