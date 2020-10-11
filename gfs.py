import os, sys, subprocess, glob, h5py, time, datetime
import pygrib
from datetime import date
import numpy as np

def readgfs(gfs_filename):
    grbs=pygrib.open(gfs_filename)
    grb = grbs.select()[:]
    
    # index별 data를 저장할 dictionary 변수 정의 
    dic_data = {}
    
    index = int(arr_var[0])
    lats,lons=grb[index].latlons() 
    
    # get the mask data according to lats, lons 
    m_po,_=lats.shape
    m_po=int(m_po/2)
    mask_lon=(lons[m_po,:] >=min_lon)&(lons[m_po,:]<=max_lon)
    mask_lat=(lats[:,m_po]>=min_lat)&(lats[:,m_po]<=max_lat)
    
    m_lon=lons[mask_lat,:]
    m_lon=m_lon[:,mask_lon]
    m_lat=lats[mask_lat,:]
    m_lat=m_lat[:,mask_lon]
    
    
        
    for i_var in arr_var: 
        #print(i_var)
        
        data_org=grb[int(i_var)].values 
        data=data_org #*3600
        #print(np.min(data), np.max(data))
  
        m_data=data[mask_lat,:]
        m_data=m_data[:,mask_lon]
        dic_data[str(i_var)] = m_data
    
        #print(m_lon.shape, m_lat.shape, m_data.shape)
        #print(np.min(data), np.max(data))
    #print(np.min(m_lon), np.max(m_lon))
    #print(dic_data)
    return m_lon, m_lat, dic_data


def w_pre_processing(out_f_name, w_lon, w_lat, w_data):
    dic_data_out = {}
    dic_data_out = w_data
    out_f_name = out_f_name.strip()
    #print(out_f_path, out_f_name)
   
    if os.path.exists(out_f_path+out_f_name+'.hdf5'):
        os.remove(out_f_path+out_f_name+'.hdf5')
    """ 
    d_lon = np.random.random(size=(i_lon, i_lat))
    d_lat = np.random.random(size=(i_lon, i_lat))
    d_data = np.random.random(size=(i_lon, i_lat))
    """
    
    f_w = h5py.File(out_f_path+out_f_name+'-gfs.hdf5', 'w')
#    f_w = h5py.File(out_f_path+'2019-gfs.hdf5', 'a')
    g_w = f_w.create_group('gfs')
    
#    g_w.create_dataset('time', data=out_f_name)
    g_w.create_dataset('lon', data=w_lon)
    g_w.create_dataset('lat', data=w_lat)
    
    for k, v in dic_data_out.items(): 
        #print(k,v)
        g_w.create_dataset(k.strip(), data=v) 
        
    f_w.close()


if __name__ == "__main__": 
    
    conf_file = open('gfs.conf', 'r')

    lines = conf_file.read()
    conf_file.close()

    lines = lines.split('\n')

    i_path = lines[0].split(':')[1].strip()
    s_date = lines[1].split(':')[1].strip()
    e_date = lines[2].split(':')[1].strip()
    i_time = int(lines[3].split(':')[1].strip())
    arr_var = lines[4].split(':')[1].strip().split(',')
    min_lon = int(lines[5].split(':')[1].strip().split(',')[0].strip())
    max_lon = int(lines[5].split(':')[1].strip().split(',')[1].strip())
    min_lat = int(lines[6].split(':')[1].strip().split(',')[0].strip())
    max_lat = int(lines[6].split(':')[1].strip().split(',')[1].strip())
    st_name = lines[7].split(':')[1].strip()
    out_f_path = lines[8].split(':')[1].strip()

    #print(i_path, s_date, e_date, min_lon, max_lon, min_lat, max_lat, st_name, out_f_path)

    s_year = date.fromisoformat(s_date).timetuple()[0]
    path = i_path + str(s_year) +'/'

    s_dates = date.fromisoformat(s_date)
    e_dates = date.fromisoformat(e_date)

    #print(s_dates, e_dates)
    #print(arr_var)
    
    s_year = date.fromisoformat(s_date).timetuple()[0]
    s_dates = date.fromisoformat(s_date)
    e_dates = date.fromisoformat(e_date)

    days = [s_dates + datetime.timedelta(days=x) for x in range((e_dates-s_dates).days +1)]
    #print(days)

    start_time = time.time() 

    for day in days: 
        #print(day.strftime('%Y%m%d'))
        s_month = day.strftime('%m')
        for i in list(range(0, 1440, i_time)): 
            h_time = int(i/60)
            m_time = int(i%60)       
            h_file = st_name + '.' + day.strftime('%Y%m%d') + str(h_time).zfill(2).strip() +'.f000'
            #print(h_file)
            path = i_path + str(s_year) +'/' + s_month +'/'
            #print(path)
            for file_name in os.listdir(path):
                if os.path.isfile(path + file_name) and h_file in file_name: 
                    #print(file_name)
                    m_lon, m_lat, m_data = readgfs(path+file_name)
                    #print(m_lon, m_lat, m_data)
                    w_pre_processing(day.strftime('%Y%m%d')+str(h_time).zfill(2)+str(m_time).zfill(2),m_lon, m_lat, m_data)

    end_time = time.time()
    #print('Working Time: {} sec'.format(end_time-start_time))                

    out_f = open(out_f_path+"result.txt", "a")
    out_f.write('gfs, Date: '+ str(time.ctime()) + ', Time Interval:' + str(i_time) +', Start Date: '+str(s_date)+', End Date: '+str(e_date)+', working time: {} sec \n'.format(end_time-start_time))
    out_f.close()

