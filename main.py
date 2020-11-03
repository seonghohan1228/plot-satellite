import os
os.environ['PROJ_LIB'] = r"C:\Users\seonghohan\anaconda3\share\proj"

from plot import plot_PC1, plot_POS, plot_MAG, plot_DET, plot_TEL

# File path is created in read_hdf.
filename_HEPD = 'HEPD_DIV_20200717_0850_ORB_08795.h5'
filename_MEPD = 'MEPD_SCI_20200717_0851_ORB_08795.h5'

# PC1 plots require both HEPD and MEPD data.
plot_PC1(filename_HEPD, filename_MEPD)

# Satellite position plot can be obtained from either HEPD or MEPD data (might vary slightly).
plot_POS(filename_HEPD)

# Magnetic field plot can be obtained from either HEPD or MEPD data (HEPD seems to be closer).
plot_MAG(filename_HEPD)

# Detector plot uses MEPD data.
plot_DET(filename_MEPD)

# Telescope plot uses HEPD data.
plot_TEL(filename_HEPD)
