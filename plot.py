import matplotlib
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
from mpl_toolkits.basemap import Basemap
from datetime import datetime
from read import read_TIME, read_DT, read_PC1, read_POS, read_MAG, read_DET, read_TEL
from calc import sliceAB, sliceAB2, mlat, B_avg

# Custom colormap, cmap
jet = plt.cm.get_cmap('jet', 256)
newcolors = jet(np.linspace(0, 1, 256))
white = np.array([256/256, 256/256, 256/256, 1])
newcolors[0, :] = white
new_cmap = matplotlib.colors.ListedColormap(newcolors)

# COMPLETE
# Plot TIME vs. PC1
def plot_PC1(PC1_HEPD, PC1_MEPD, TIME_HEPD, TIME_MEPD, DT):
    # Divide PC1 data into MEPD-A and MEPD-B
    PC1_MEPD_A, PC1_MEPD_B = sliceAB(PC1_MEPD, DT, 3)
    TIME_MEPD_A, TIME_MEPD_B = sliceAB(TIME_MEPD, DT, 3)

    # Plot PC1 data.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Plot HEPD PC1 data.
    ax1.plot(PC1_HEPD, TIME_HEPD - TIME_HEPD[0], '-k')

    # Plot MEPD-A and MEPD-B PC1 data.
    ax2.plot(PC1_MEPD_A, TIME_MEPD_A - TIME_MEPD_A[0], '-k', label='MEPD-A')
    ax2.plot(PC1_MEPD_B, TIME_MEPD_B - TIME_MEPD_B[0], '-r', label='MEPD-B')

    ax1.set_title('HEPD: Time vs PC1')
    ax1.set_xlabel('PC1')
    ax1.set_ylabel('Time (sec)')

    ax2.set_title('MEPD: Time vs PC1')
    ax2.set_xlabel('PC1')

    plt.legend()

    plt.savefig('./plots/PC1.png')
    plt.show()
    plt.close('all')


# INCOMPLETE: Geomagnetic lattitude plot
# Plot satellite position onto map
def plot_POS(POS, TIME, POLE):
    LAT = POS[:, 0]
    LON = POS[:, 1]
    #ALT = POS[:, 2]
    
    start_time = datetime.fromtimestamp(TIME[0]) # Converts UNIX datetime to UTC datetime
    end_time = datetime.fromtimestamp(TIME[-1])

    # Draw maps.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))

    # Mercador projection.
    m1 = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85, llcrnrlon=-180,urcrnrlon=180, ax=ax1)
    m1.drawcoastlines()
    m1.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    m1.drawparallels(np.arange(-90,90,30), labels=[False, True, False, False])

    # Orthographic projection.
    m2 = Basemap(projection='ortho', lat_0=POLE, lon_0=0, ax=ax2)
    m2.drawcoastlines()
    m2.drawmeridians(np.arange(0,360,45), labels=[False, False, False, True])
    m2.drawparallels(np.arange(-90,90,30)) # Cannot label parallels on Orthographic basemap.
    
    # Plot satellite position.
    X, Y = m1(LON, LAT)
    s1 = m1.scatter(X, Y, marker='.', c=TIME, cmap=plt.cm.jet)
    ax1.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax1.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))

    X, Y = m2(LON, LAT)
    s2 = m2.scatter(X, Y, marker='.', c=TIME, cmap=plt.cm.jet)
    ax2.annotate(start_time.strftime('%H:%M'), (X[0], Y[0]))
    ax2.annotate(end_time.strftime('%H:%M'), (X[-1], Y[-1]))
    
    cbar = fig.colorbar(s1, ax=ax2, label='Time (UNIX)')

    # Plot terminator
    m1.nightshade(start_time)
    m2.nightshade(end_time)

    ax1.set_title('Orbit (Mercador projection)')
    ax2.set_title('Orbit (Orthographic projection)')

    plt.savefig('./plots/Position.png')
    plt.show()
    plt.close('all')
    

# COMPLETE
# Plot Magnetic Field vs. TIME.
def plot_MAG(filename):
    # Read TIME data
    TIME = read_TIME(filename)

    # Read magnetic field data
    MAG = read_MAG(filename)
    MAG_avg = B_avg(MAG)

    fig, ax = plt.subplots()
    ax.plot(TIME, MAG[:, 0], 'k', label='Bx')
    ax.plot(TIME, MAG[:, 1], 'b', label='By')
    ax.plot(TIME, MAG[:, 2], 'r', label='Bz')
    ax.plot(TIME, MAG[:, 4], '--k', label='IGRF Bx')
    ax.plot(TIME, MAG[:, 5], '--b', label='IGRF By')
    ax.plot(TIME, MAG[:, 6], '--r', label='IGRF Bz')
    ax.plot(TIME, MAG_avg, '--y', label='IGRF|B|')
    
    ax.tick_params(which='both', direction='in')
    #plt.xlabel('Time')
    plt.ylabel('Magnetic Field (nT)')
    plt.ylim(-60000, 60000)
    plt.legend(loc='upper center', ncol=7, prop={'size': 6})

    plt.savefig('./plots/Magnetic Field.png')
    #plt.show()
    plt.close('all')


# INCOMPLETE: y-axis, colorbar
# Plot Detector Data (Energy) vs. TIME
def plot_DET(filename):
    # Read TIME data
    TIME = read_TIME(filename)

    # Read detector data
    DET0, DET1, DET2, DET3 = read_DET(filename)
    
    # Read DT (Subunit ID) data
    DT = read_DT(filename)
    
    # Divide PC1 data into MEPD-A and MEPD-B
    DET0A, DET0B = sliceAB2(DET0, DT, 3)
    DET1A, DET1B = sliceAB2(DET1, DT, 3)
    DET2A, DET2B = sliceAB2(DET2, DT, 3)
    DET3A, DET3B = sliceAB2(DET3, DT, 3)
    TIMEA, TIMEB = sliceAB(TIME, DT, 3)
    
    # Plot MEPD-A
    xmin = mdates.date2num(datetime.fromtimestamp(TIMEA[0]))
    xmax = mdates.date2num(datetime.fromtimestamp(TIMEA[-1]))
    
    figA, (axA0, axA1, axA2, axA3) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(8, 6))
    imA = axA0.imshow(
        X=np.transpose(DET0A), 
        origin='lower', cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imA = axA1.imshow(
        X=np.transpose(DET1A), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imA = axA2.imshow(
        X=np.transpose(DET2A), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imA = axA3.imshow(
        X=np.transpose(DET3A), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])

    axA0.text(0.02, 0.85, 'Detector 0', horizontalalignment='left', transform=axA0.transAxes)
    axA1.text(0.02, 0.85, 'Detector 1', horizontalalignment='left', transform=axA1.transAxes)
    axA2.text(0.02, 0.85, 'Detector 2', horizontalalignment='left', transform=axA2.transAxes)
    axA3.text(0.02, 0.85, 'Detector 3', horizontalalignment='left', transform=axA3.transAxes)
    
    figA.subplots_adjust(left=0.08, right=0.9, bottom=0.05, top=0.95, wspace=0.0, hspace=0.0)
    cb_axA = figA.add_axes([0.92, 0.05, 0.02, 0.9])
    cbarA = figA.colorbar(imA, cax=cb_axA)
    axA1.set_ylabel('Energy [keV]') # Divde y label into two parts for center alignment
    axA2.set_ylabel('MEPD-A')
    axA3.xaxis_date()
    axA3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    plt.savefig('./plots/MEPD-A Detector.png')
    
    # Plot MEPD-B
    figB, (axB0, axB1, axB2, axB3) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(8, 6))
    imB = axB0.imshow(
        X=np.transpose(DET0B), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imB = axB1.imshow(
        X=np.transpose(DET1B), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imB = axB2.imshow(
        X=np.transpose(DET2B), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])
    imB = axB3.imshow(
        X=np.transpose(DET3B), 
        origin='lower', 
        cmap=new_cmap, 
        aspect='auto', 
        interpolation='none', 
        extent = [xmin, xmax, 0, 400])

    axB0.text(0.02, 0.85, 'Detector 0', horizontalalignment='left', transform=axB0.transAxes)
    axB1.text(0.02, 0.85, 'Detector 1', horizontalalignment='left', transform=axB1.transAxes)
    axB2.text(0.02, 0.85, 'Detector 2', horizontalalignment='left', transform=axB2.transAxes)
    axB3.text(0.02, 0.85, 'Detector 3', horizontalalignment='left', transform=axB3.transAxes)

    figB.subplots_adjust(left=0.08, right=0.9, bottom=0.05, top=0.95, wspace=0.0, hspace=0)
    cb_axB = figB.add_axes([0.92, 0.05, 0.02, 0.9])
    cbarB = figB.colorbar(imB, cax=cb_axB)
    axB1.set_ylabel('Energy [keV]')
    axB2.set_ylabel('MEPD-B')
    axB3.xaxis_date()
    axB3.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))

    plt.savefig('./plots/MEPD-B Detector.png')
    #plt.show()
    plt.close('all')

# INCOMPLETE: HEPD proton and electron not divided
# Plot Telescope Data (Energy) vs. TIME
def plot_TEL(filename):
    # Read TIME data
    TIME = read_TIME(filename)

    # Read detector data
    TEL0, TEL1, TEL2 = read_TEL(filename)
    
    # Plot MEPD-A
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 6))
    im = ax0.imshow(X=np.transpose(TEL0), origin='lower', cmap=new_cmap, aspect='auto', interpolation='none')
    im = ax1.imshow(X=np.transpose(TEL1), origin='lower', cmap=new_cmap, aspect='auto', interpolation='none')
    im = ax2.imshow(X=np.transpose(TEL2), origin='lower', cmap=new_cmap, aspect='auto', interpolation='none')
    ax0.text(0.02, 0.85, 'Telescope 0', horizontalalignment='left', transform=ax0.transAxes)
    ax1.text(0.02, 0.85, 'Telescope 1', horizontalalignment='left', transform=ax1.transAxes)
    ax2.text(0.02, 0.85, 'Telescope 2', horizontalalignment='left', transform=ax2.transAxes)
    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.05, top=0.95, wspace=0.0, hspace=0.0)
    cb_ax = fig.add_axes([0.92, 0.05, 0.02, 0.9])
    cbar = fig.colorbar(im, cax=cb_ax)
    ax1.set_ylabel('Energy [keV]')

    plt.savefig('./plots/HEPD Telescope.png')
    #plt.show()
    plt.close('all')
