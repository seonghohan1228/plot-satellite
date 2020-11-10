# Required for proj_lib to work properly with basemap.
# Change path the the package location of proj.
import os
os.environ['PROJ_LIB'] = r"C:\Users\seonghohan\anaconda3\share\proj"

from read import read_hdf
from select import select_HEPD, select_MEPD
from plot import plot_PC1, plot_POS, plot_MAG, plot_DET, plot_TEL

# Filename (file path is created automatically in hdf_read()).
HEPD_FILENAME = 'HEPD_DIV_20200717_0850_ORB_08795.h5'
MEPD_FILENAME = 'MEPD_SCI_20200717_0851_ORB_08795.h5'

NORTH_POLE = 90
SOUTH_POLE = -90


# Read and store data.
HEPD_dataset1, HEPD_dataset2 = read_hdf(HEPD_FILENAME)
MEPD_dataset1, MEPD_dataset2 = read_hdf(MEPD_FILENAME)

# Select required data from datasets.
HEPD_time, HEPD_pc1, HEPD_pos, HEPD_mag, tel0, tel1, tel2 = select_HEPD(HEPD_dataset1, HEPD_dataset2)
MEPD_time, dt, MEPD_pc1, MEPD_pos, MEPD_mag, det0, det1, det2, det3 = select_MEPD(MEPD_dataset1, MEPD_dataset2)

# PC1 plots require both HEPD and MEPD data.
#plot_PC1(HEPD_pc1, MEPD_pc1, HEPD_time, MEPD_time, dt)

# Satellite position plot can be obtained from either HEPD or MEPD data (might vary slightly).
plot_POS(MEPD_pos, MEPD_time, SOUTH_POLE)

# Magnetic field plot can be obtained from either HEPD or MEPD data (HEPD seems to be closer).
plot_MAG(HEPD_FILENAME)

# Detector plot uses MEPD data.
plot_DET(MEPD_FILENAME)

# Telescope plot uses HEPD data.
plot_TEL(HEPD_FILENAME)
