import numpy as np
import pyIGRF


# Slice data into MEPD-A and MEPD-B (works for PC1 and TIME)
def sliceAB(data, subunit, n):
    A = []
    B = []
    for i in range(len(data)):
        if subunit[i] == n:
            A.append(data[i])
        else:
            B.append(data[i])

    return A, B


# sliceAB 2D array data of MEPD data
def sliceAB2(data, subunit, n):
    A = []
    B = []
    for i in range(len(data)):
        if subunit[i] == n:
            A.append(data[i, :])
        else:
            B.append(data[i, :])
    
    return A, B

"""
# Calculate pole from IGRF data
def pole(year, alt):
    Pole = np.zeros((181, 360))
    Nmin, Smin = 5, -5
    Nlat, Nlon, Slat, Slon = 0, 0, 0, 0

    for i in range(180):
        for j in range(91):
            """"""
            declination (+ve east)
            inclination (+ve down)
            horizontal intensity
            north, east component
            vertical component (+ve down)
            total intensity unit: degreee or nT
            """"""
            D, I, H, X, Y, X, F = pyIGRF.igrf_value(2*j - 90, 2*i, alt, year)
            if I > 0 and I < Nmin:
                Nmin = I
                Nlat, Nlon = j, i
            if I < 0 and I > Smin:
                Smin = I
                Slat, Slon = j, i
    
    return Nlat, Nlon, Slat, Slon
"""

# Calculate geomagnetic latitude data from pole
def mlat(N_mag_pole, S_mag_pole):
    pass

# Calculate average magnetic field.
def B_avg(B):
    avg = []
    for i in range(len(B)):
        avg.append(np.sqrt((B[i, 4])**2 + (B[i, 5])**2 + (B[i, 6])**2))

    return avg
