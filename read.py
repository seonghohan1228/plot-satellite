import h5py
import numpy as np


# Read TIME data
def read_TIME(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    if group == 'HEPD_DIV':
        data = dataset[:, 7]
    if group == 'MEPD_SCI':
        data = dataset[:, 10]

    return data


# Read DT (Subunit ID) data
def read_DT(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data = dataset[:, 4]

    return data


# Read PC1 data
def read_PC1(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data = dataset[:, 5]

    return data


# Read satellite position data
def read_POS(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block2_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data = dataset[:, 16:19]

    return data


# Read magnetic field data
def read_MAG(filename):
    filepath = 'data/' + filename

    # group: HEPD_DIV or MEPD_SCI
    group = filename[0:8]
    dataset = 'block2_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data = dataset[:, 0:8]

    return data


# Read detector data
def read_DET(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data0 = dataset[:, 13:77]
    data1 = dataset[:, 81:145]
    data2 = dataset[:, 149:213]
    data3 = dataset[:, 217:281]

    return data0, data1, data2, data3


# Read telescope data
def read_TEL(filename):
    filepath = 'data/' + filename

    # group: HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = np.array(hdf[path])
    
    data0 = dataset[:, 9:49] # Last values are checksums
    data1 = dataset[:, 50:90]
    data2 = dataset[:, 91:131]

    return data0, data1, data2
