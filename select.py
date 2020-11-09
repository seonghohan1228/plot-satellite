## select.py ##
## Functions that select required data from given datasets.

# Data index
HEPD_TIME = 7       # Time (UNIX)
MEPD_TIME = 10      # Time (UNIX)
DT = 4              # Subunit ID
PC1 = 5             # Packet count 1
POS = 16            # Position (deg)
POS_LEN = 3         # Position data length
MAG = 0             # Magnetic field
MAG_LEN = 8         # Magnetic field data length
DET0 = 13           # Detector 0
DET1 = 81           # Detector 1
DET2 = 149          # Detector 2
DET3 = 217          # Detector 3
DET_LEN = 64        # Detector data length
TEL0 = 9            # Telescope 0
TEL1 = 50           # Telescope 1
TEL2 = 91           # Telescope 2
TEL_LEN = 40        # Telescope data length


# Functions
def select_HEPD(dataset1, dataset2):
    time = dataset1[:, HEPD_TIME]
    pc1 = dataset1[:, PC1]
    pos = dataset2[:, POS:POS + POS_LEN]
    mag = dataset2[:, MAG:MAG + MAG_LEN]
    tel0 = dataset1[:, TEL0:TEL0 + TEL_LEN]
    tel1 = dataset1[:, TEL1:TEL1 + TEL_LEN]
    tel2 = dataset1[:, TEL2:TEL2 + TEL_LEN]

    return time, pc1, pos, mag, tel0, tel1, tel2


def select_MEPD(dataset1, dataset2):
    time = dataset1[:, MEPD_TIME]
    dt = dataset1[:, DT]
    pc1 = dataset1[:, PC1]
    pos = dataset2[:, POS:POS + POS_LEN]
    mag = dataset2[:, MAG:MAG + MAG_LEN]
    det0 = dataset1[:, DET0:DET0 + DET_LEN]
    det1 = dataset1[:, DET1:DET1 + DET_LEN]
    det2 = dataset1[:, DET2:DET2 + DET_LEN]
    det3 = dataset1[:, DET3:DET3 + DET_LEN]

    return time, dt, pc1, pos, mag, det0, det1, det2, det3

