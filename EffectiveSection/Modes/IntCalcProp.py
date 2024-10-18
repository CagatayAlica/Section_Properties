import math
import numpy as np


def calcProps(data: []):
    """
    This function calculates the section properties using effective length and thicknesses.
    :param data: 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
    :return: Ar, ygc, Ix, xgc, Iy
    """
    # data:
    # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
    # Area of cross-section parts
    array_Ar = []
    for i in data:
        dai = i[5] * math.sqrt(math.pow(i[4] - i[2], 2) + math.pow(i[3] - i[1], 2))
        array_Ar.append(dai)
    array_Ar = np.array(array_Ar)
    # =======
    Ar = np.sum(array_Ar)
    # =======
    # First moment of area X
    array_Sx0 = []
    for i in data:
        sx0i = (i[2] + i[4]) * (i[5] * math.sqrt(math.pow(i[4] - i[2], 2) + math.pow(i[3] - i[1], 2))) / 2
        array_Sx0.append(sx0i)
    array_Sx0 = np.array(array_Sx0)
    Sx0 = np.sum(array_Sx0)
    # =======
    ygc = Sx0 / Ar
    # =======
    # Second moment of area X
    array_Ix0 = []
    for i in data:
        Ix0 = (math.pow(i[2], 2) + math.pow(i[4], 2) + i[2] * i[4]) * (
                i[5] * math.sqrt(math.pow(i[4] - i[2], 2) + math.pow(i[3] - i[1], 2)) / 3)
        array_Ix0.append(Ix0)
    array_Ix0 = np.array(array_Ix0)
    Ix0 = np.sum(array_Ix0)
    # =======
    Ix = Ix0 - Ar * math.pow(ygc, 2)
    # =======

    # First moment of area Y
    array_Sy0 = []
    for i in data:
        sx0i = (i[1] + i[3]) * (i[5] * math.sqrt(math.pow(i[4] - i[2], 2) + math.pow(i[3] - i[1], 2))) / 2
        array_Sy0.append(sx0i)
    array_Sy0 = np.array(array_Sy0)
    Sy0 = np.sum(array_Sy0)
    # =======
    xgc = Sy0 / Ar
    # =======
    # Second moment of area Y
    array_Iy0 = []
    for i in data:
        Iy0 = (math.pow(i[1], 2) + math.pow(i[3], 2) + i[1] * i[3]) * (
                i[5] * math.sqrt(math.pow(i[4] - i[2], 2) + math.pow(i[3] - i[1], 2)) / 3)
        array_Iy0.append(Iy0)
    array_Iy0 = np.array(array_Iy0)
    Iy0 = np.sum(array_Iy0)
    # =======
    Iy = Iy0 - Ar * math.pow(xgc, 2)
    # =======
    return Ar, ygc, Ix, xgc, Iy
