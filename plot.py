import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from read import *
from calc import *

# UNFINISHED: cannot differentiate MEPD-A and MEPD-B, missing plot details
# Plot TIME vs. PC1
def plot_PC1(filename):
    # Read TIME data
    TIME = read_TIME(filename)
    # Read PC1 data
    PC1 = read_PC1(filename)

    fig, ax = plt.subplots()
    ax.plot(PC1, TIME, '.')
    plt.show()


# UNFINISHED: cannot plot onto map
# Plot satellite position onto map
def plot_POS(filename):
    # Read POS data
    POS = read_POS(filename)
    LAT = POS[:, 0]
    LONG = POS[:, 1]
    


# UNFINISHED: missing plot details
# Plot Magnetic Field vs. TIME.
def plot_MAG(filename):
    # Read TIME data
    TIME = read_TIME(filename)

    # Read magnetic field data
    MAG = read_MAG(filename)
    MAG_avg = B_avg(MAG)

    fig, ax = plt.subplots()
    for i in range(8):
        if i < 3:
            ax.plot(TIME, MAG[:, i])
        if 3 < i < 7:
            ax.plot(TIME, MAG[:, i], '--')
    ax.plot(TIME, MAG_avg, '--')
    plt.show()







