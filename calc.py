import numpy


# Slice data into MEPD-A and MEPD-B (works for PC1 and TIME)
def slicePC1(data):
    A = []
    B = []
    for i in range(len(data)):
        if i % 2:
            A.append(data[i])
        else:
            B.append(data[i])

    return A, B


# Calculate average magnetic field.
def B_avg(B):
    avg = []
    for i in range(len(B)):
        avg.append(numpy.sqrt((B[i, 4])**2 + (B[i, 5])**2 + (B[i, 6])**2))

    return avg
