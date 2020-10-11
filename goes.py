import os, sys, subprocess, glob, h5py, time, datetime, sys
from datetime import date
from netCDF4 import Dataset
import numpy as np
from ipykernel import kernelapp as app
import tracemalloc
import time as t 
import osgeo import osr
import osgeo import gdal

# Install these libraries before exec python file 
# goes.conf 과 goes.py 파일은 같은 디렉터리에 있어야 함 
# conda install theano keras graphviz numpy scipy scikit-learn matplotlib pandas pydot h5py gdal -y
# python goes.py 2019-01-01 2019-12-31 


# Define KM_PER_DEGREE
KM_PER_DEGREE = 111.32

# GOES-16 Spatial Reference System
sourcePrj = osr.SpatialReference()
sourcePrj.ImportFromProj4('+proj=geos +h=35786023.0 +a=6378137.0 +b=6356752.31414 +f=0.0033528106811935607 +lat_0=0.0 +lon_0=-137 +sweep=x +no_defs')

# Lat/lon WSG84 Spatial Reference System
targetPrj = osr.SpatialReference()
targetPrj.ImportFromProj4('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs')

def exportImage(image,path):
    driver = gdal.GetDriverByName('netCDF')
    return driver.CreateCopy(path,image,0)

def getGeoT(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]

def getScaleOffset(path):
    nc = Dataset(path, mode='r')
    scale = nc.variables['Rad'].scale_factor
    offset = nc.variables['Rad'].add_offset
    nc.close()
    return scale, offset
    
def remap(path, extent, resolution, x1, y1, x2, y2):
    
    # GOES-16 Extent (satellite projection) [llx, lly, urx, ury]
    GOES16_EXTENT = [x1, y1, x2, y2]
    
    # Setup NetCDF driver
    gdal.SetConfigOption('GDAL_NETCDF_BOTTOMUP', 'NO')
        
    # Read scale/offset from file
    scale, offset = getScaleOffset(path) 
         
    try:  
        connectionInfo = 'NETCDF:\"' + path + '\":Rad'    
        # Open NetCDF file (GOES-16 data)  
        raw = gdal.Open(connectionInfo)
    except:
        connectionInfo = 'HDF5:\"' + path + '\"://Rad'    
        # Open NetCDF file (GOES-16 data)  
        raw = gdal.Open(connectionInfo)    
                
    # Setup projection and geo-transformation
    raw.SetProjection(sourcePrj.ExportToWkt())
    #raw.SetGeoTransform(getGeoT(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))
    raw.SetGeoTransform(getGeoT(GOES16_EXTENT, raw.RasterYSize, raw.RasterXSize))  
  
    #print (KM_PER_DEGREE)
    # Compute grid dimension
    sizex = int(((extent[2] - extent[0]) * KM_PER_DEGREE) / resolution)
    sizey = int(((extent[3] - extent[1]) * KM_PER_DEGREE) / resolution)
    
    # Get memory driver
    memDriver = gdal.GetDriverByName('MEM')
   
    # Create grid
    grid = memDriver.Create('grid', sizex, sizey, 1, gdal.GDT_Float32)
    
    # Setup projection and geo-transformation
    grid.SetProjection(targetPrj.ExportToWkt())
    grid.SetGeoTransform(getGeoT(extent, grid.RasterYSize, grid.RasterXSize))

    # Perform the projection/resampling 

    #print ('Remapping', path)
        
    start = t.time()
   
    gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=ALL_CPUS']) 
    #gdal.ReprojectImage(raw, grid, sourcePrj.ExportToWkt(), targetPrj.ExportToWkt(), gdal.GRA_NearestNeighbour, options=['NUM_THREADS=40']) 
    
    #print ('- finished! Time:', t.time() - start, 'seconds')
    
    # Close file
    raw = None
        
    # Read grid data
    array = grid.ReadAsArray()
    
    # Mask fill values (i.e. invalid values)
    np.ma.masked_where(array, array == -1, False)
    
    # Apply scale and offset
    array = array * scale + offset
    
    grid.GetRasterBand(1).SetNoDataValue(-1)
    grid.GetRasterBand(1).WriteArray(array)

    return grid


def exportImage(image,path):
    driver = gdal.GetDriverByName('netCDF')
    return driver.CreateCopy(path,image,0)

def getGeoT(extent, nlines, ncols):
    # Compute resolution based on data dimension
    resx = (extent[2] - extent[0]) / ncols
    resy = (extent[3] - extent[1]) / nlines
    return [extent[0], resx, 0, extent[3] , 0, -resy]

def getScaleOffset(path):
    nc = Dataset(path, mode='r')
    scale = nc.variables['Rad'].scale_factor
    offset = nc.variables['Rad'].add_offset
    nc.close()
    return scale, offset
    

def readgoes(goes_filename, ch):
    #print(goes_filename)
    g17nc=Dataset(goes_filename,'r')
    #radiance=g17nc.variables[attrs[0].strip()][:]
    pla_fk1=g17nc.variables[attrs[1].strip()][:]
    pla_fk2=g17nc.variables[attrs[2].strip()][:]
    pla_bc1=g17nc.variables[attrs[3].strip()][:]
    pla_bc2=g17nc.variables[attrs[4].strip()][:]
    kappa=g17nc.variables[attrs[5].strip()][:]

    extent = [min_lon, min_lat, max_lon, max_lat]
 
    # Choose the image resolution (the higher the number the faster the processing is)
    resolution = 2
 
    # Calculate the image extent required for the reprojection
    H = g17nc.variables['goes_imager_projection'].perspective_point_height
    x1 = g17nc.variables['x_image_bounds'][0] * H
    x2 = g17nc.variables['x_image_bounds'][1] * H
    y1 = g17nc.variables['y_image_bounds'][1] * H
    y2 = g17nc.variables['y_image_bounds'][0] * H
 
    # Call the reprojection funcion
    grid = remap(goes_filename, extent, resolution, x1,y1,x2,y2)
    radiance = grid.ReadAsArray()

    #radiance to data
    if ch <= 6:
        mask_data=kappa*radiance
    else:
        mask_data=(pla_fk2/(np.log(pla_fk1/radiance+1))-pla_bc1)/pla_bc2
    #print(str(ch), data.shape, flush=True)
    #print(np.min(data), np.max(data),flush=True)
    #print(np.min(radiance), np.max(radiance),flush=True)

    grid_x, grid_y = np.mgrid[-125:-105:1113j,30:10:1113j]
    #grid_x=np.flip(grid_x,axis=1)
    grid_y=np.flip(grid_y,axis=1)
    #print(grid_x.shape, grid_y.shape)
    grid_x=grid_x+360
    #print(np.min(grid_x), np.max(grid_x))
    mask_data=np.flip(mask_data,axis=0)

    return mask_data.T, grid_x, grid_y    

def w_pre_processing(out_f_name, w_lon, w_lat, w_data):
    out_f_name = out_f_name.strip()
    dic_data_out = {} 
    dic_data_out = w_data
    
    if os.path.exists(out_f_path+out_f_name+'-goes.hdf5'):
        os.remove(out_f_path+out_f_name+'-goes.hdf5')
    
    d_lon = np.random.random(size=(i_lon, i_lat))
    d_lat = np.random.random(size=(i_lon, i_lat))
    d_data = np.random.random(size=(i_lon, i_lat))
    
    f_w = h5py.File(out_f_path+out_f_name+'-goes.hdf5', 'w')
    g_w = f_w.create_group('goes')
    
    g_w.create_dataset('lon', data=w_lon)
    g_w.create_dataset('lat', data=w_lat)
    for k, v in dic_data_out.items():
        g_w.create_dataset(k.strip(), data=v)
    
    f_w.close()

if __name__ == "__main__": 
    #os.environ["OMP_NUM_THREADS"] = "10"
    # read goes.conf file
    conf_file = open('goes.conf', 'r')

    lines = conf_file.read()
    conf_file.close()

    lines = lines.split('\n')

    i_path = lines[0].split(':')[1].strip()
    #s_date = lines[1].split(':')[1].strip()
    #e_date = lines[2].split(':')[1].strip()
    
    # for parallel python exec input pars 
    s_date = sys.argv[1].strip()
    e_date = sys.argv[2].strip()
    i_time = int(lines[3].split(':')[1].strip())
    i_lon = 5423
    i_lat = 5423
    min_lon = int(lines[4].split(':')[1].strip().split(',')[0].strip())
    max_lon = int(lines[4].split(':')[1].strip().split(',')[1].strip())
    min_lat = int(lines[5].split(':')[1].strip().split(',')[0].strip())
    max_lat = int(lines[5].split(':')[1].strip().split(',')[1].strip())
    attrs = lines[6].split(':')[1].strip().split(',')

    st_name = lines[7].split(':')[1].strip()
    out_f_path = lines[8].split(':')[1].strip()
    channels = lines[9].split(':')[1].strip().split(',')
    #channel = int(sys.argv[1].strip())
    #print(channel)
    s_year = date.fromisoformat(s_date).timetuple()[0]
    s_days = int(date.fromisoformat(s_date).timetuple()[7])
    e_days = int(date.fromisoformat(e_date).timetuple()[7])
    
    if max_lon > 180:
        #mlon=max_lon
        max_lon=max_lon-360

    if min_lon > 180:
        min_lon=min_lon-360
    

    # print(channel)
    start_c_time = (time.ctime())
    start_time = time.time() 
    tracemalloc.start() 

    s_dates = date.fromisoformat(s_date)
    e_dates = date.fromisoformat(e_date)

    days = [s_dates + datetime.timedelta(days=x) for x in range((e_dates-s_dates).days +1)]
    # print(days)

    for day in days: 
        for i in list(range(0, 1440, i_time)):         
            h_time = int(i/60)
            m_time = int(i%60)
            path = i_path + str(s_year) + '/' + str(date.fromisoformat(str(day)).timetuple()[7]).zfill(3) + '/' + str(h_time).zfill(2) +'/'
            dic_data = {} 
            m_lon = []
            m_lat = [] 
            for c in channels: 
                h_file = c.zfill(2).strip() + '_G17' + '_' +'s'+str(s_year)+str(date.fromisoformat(str(day)).timetuple()[7]).zfill(3)+str(h_time).zfill(2)+str(m_time).zfill(2)
                if not os.path.exists(path): 
                    #print('dir %s does not exists'%path)
                    os.mkdir(path)
                else: 
                    for file_name in os.listdir(path): 
                        if os.path.isfile(path+file_name) and h_file in file_name:
                            #print(file_name)
                            
                            dic_data[c], m_lon, m_lat = readgoes(path+file_name, int(c))

            w_pre_processing(day.strftime('%Y%m%d') + str(h_time).zfill(2) + str(m_time).zfill(2), m_lon, m_lat, dic_data )
            

    end_time = time.time()
    #print('Working Time: {} sec'.format(end_time-start_time))
    current, peak = tracemalloc.get_traced_memory()
    out_f = open(out_f_path+"result.txt", "a")
    out_f.write('Start: ' + start_c_time + ', End: '+ str(time.ctime())+ ', Time Interval:' + str(i_time) +', Start Date: '+str(s_date)+', End Date: '+str(e_date)+', working time: {} sec'.format(end_time-start_time) + f", Current memory usage is {current / 10**6}MB, Peak was {peak / 10**6}MB\n")
    tracemalloc.stop()
    out_f.close()

