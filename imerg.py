import os, sys, subprocess, glob, time
from datetime import date
import datetime, h5py
import numpy as np
import pandas
import tracemalloc

# 데이터 전처리 과정 파일 읽기 
def pre_processing(p_path):
    #print(p_path)
    f = h5py.File(p_path, 'r')
    group = f['Grid']

    """
    for key in group.keys():
        print('')
    """    
    attr_lon = attrs[0].strip()
    attr_lat = attrs[1].strip()
    attr_pre = attrs[2].strip() 
    
    n_lon = np.array(group[attr_lon][:])
    n_lat = np.array(group[attr_lat][:])
    n_precipitationCal = group[attr_pre][:]
    precip = np.transpose(n_precipitationCal)


    x, y = np.float32(np.meshgrid(n_lon, n_lat))
    mask_lon = (x[i_lon, :] >= min_lon)&(x[i_lon, :]<=max_lon)
    mask_lat = (y[:, i_lat] >= min_lat)&(y[:,i_lat]<=max_lat)
    
    mask_data = precip[mask_lat, :]
    mask_data = mask_data[:, mask_lon]
    m_data = np.squeeze(mask_data, axis=2)
    
    m_lon = x[mask_lat, :]
    m_lon = m_lon[:, mask_lon]+360
    m_lat = y[mask_lat,]
    m_lat = m_lat[:, mask_lon] 
    """
    print(m_data.shape, m_lon.shape, m_lat.shape)
    print(np.min(m_data), np.max(m_data))
    print(mask_data.shape)
    print(np.min(m_lon),np.max(m_lon)) #m_lat, m_data)
    """
    f.close()
    return m_lon, m_lat, m_data

# 데이터 전처리 과정 파일 쓰기 
def w_pre_processing(out_f_name, w_lon, w_lat, w_data):
    out_f_name = out_f_name.strip()
    #print(out_f_path, out_f_name)
    #print('writing data to hdf5 file') 
    if os.path.exists(out_f_path+out_f_name+'.hdf5'):
        os.remove(out_f_path+out_f_name+'.hdf5')
    
    d_lon = np.random.random(size=(i_lon, i_lat))
    d_lat = np.random.random(size=(i_lon, i_lat))
    d_data = np.random.random(size=(i_lon, i_lat))
    
    f_w = h5py.File(out_f_path+out_f_name+'-imerg.hdf5', 'w')
    g_w = f_w.create_group('imerg')
    
    g_w.create_dataset('lon', data=w_lon)
    g_w.create_dataset('lat', data=w_lat)
    g_w.create_dataset('data', data=w_data)
    
    f_w.close()

if __name__ == "__main__": 
    conf_file = open('region3_imerg.conf', 'r')

    lines = conf_file.read()
    conf_file.close()

    lines = lines.split('\n')

    i_path = lines[0].split(':')[1].strip()
    #s_date = lines[1].split(':')[1].strip()
    #e_date = lines[2].split(':')[1].strip()
    s_date = sys.argv[1].strip()
    e_date = sys.argv[2].strip()
    i_time = int(lines[3].split(':')[1].strip())
    i_lon = 900
    i_lat = 1800
    min_lon = int(lines[4].split(':')[1].strip().split(',')[0].strip())
    max_lon = int(lines[4].split(':')[1].strip().split(',')[1].strip())
    min_lat = int(lines[5].split(':')[1].strip().split(',')[0].strip())
    max_lat = int(lines[5].split(':')[1].strip().split(',')[1].strip())
    attrs = lines[6].split(':')[1].strip().split(',')
    st_name = lines[7].split(':')[1].strip()
    out_f_path = lines[8].split(':')[1].strip()

    if max_lon > 180:
        #mlon=max_lon
        max_lon=max_lon-360

    if min_lon > 180:
        min_lon=min_lon-360

        
        
    start_time = time.time() 
    tracemalloc.start() 

    s_year = date.fromisoformat(s_date).timetuple()[0]
    path = i_path + str(s_year) +'/'

    s_dates = date.fromisoformat(s_date)
    e_dates = date.fromisoformat(e_date)


    days = [s_dates + datetime.timedelta(days=x) for x in range((e_dates-s_dates).days +1)]
    for day in days: 
        # print(day.strftime('%Y%m%d'))

        for i in list(range(0, 1440, i_time)): 
            h_time = int(i/60)
            m_time = int(i%60)
            h_file = st_name + '.' + day.strftime('%Y%m%d') + '-S' + str(h_time).zfill(2) + str(m_time).zfill(2)            
            for file_name in os.listdir(path):
                if os.path.isfile(path+file_name) and h_file in file_name: 
                    #print(file_name)
                    # read data from imerg 
                    #print(path+file_name)
                    m_lon, m_lat, m_data = pre_processing(path+file_name)
                    w_pre_processing(day.strftime('%Y%m%d') + str(h_time).zfill(2)+str(m_time).zfill(2), m_lon, m_lat, m_data )
                    #print(m_lon, m_lat, m_data)

    end_time = time.time()
    #print('Working Time: {} sec'.format(end_time-start_time))
    
    current, peak = tracemalloc.get_traced_memory()
    out_f = open(out_f_path+"result.txt", "a")
    out_f.write('imerg, Date: '+ str(time.ctime()) + ', Time Interval:' + str(i_time) +', Start Date: '+str(s_date)+', End Date: '+str(e_date)+', working time: {} sec '.format(end_time-start_time) + f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB\n")
    tracemalloc.stop()
    
    out_f.close()    
