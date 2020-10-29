import h5py
import numpy


# Read TIME data
def read_TIME(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = numpy.array(hdf[path])
    
    if group == 'HEPD_DIV':
        data = dataset[:, 7]
    if group == 'MEPD_SCI':
        data = dataset[:, 10]

    return data


# Read PC1 data
def read_PC1(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block1_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = numpy.array(hdf[path])
    
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
        dataset = numpy.array(hdf[path])
    
    data = dataset[:, 16:18]

    return data


# Read magnetic field data
def read_MAG(filename):
    filepath = 'data/' + filename

    # group: MEPD_SCI or HEPD_DIV
    group = filename[0:8]
    dataset = 'block2_values'

    with h5py.File(filepath, 'r') as hdf:
        path = '/'+ group + '/' + dataset
        dataset = numpy.array(hdf[path])
    
    data = dataset[:, 0:8]

    return data


