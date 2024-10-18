import math


def lamp(b: float, t: float, ksig: float, scomed: float, iscomp: bool):
    """

    :param iscomp: is element under compression. True if compression.
    :param b: element length
    :param t: element thickness
    :param ksig: bukcling factor
    :param scomed: stress on element
    :return: lambda
    """
    if iscomp:
        lamp = (b / t) / (28.4 * math.sqrt(235.0 / scomed) * math.sqrt(ksig))
    else:
        lamp = 0.0
    return lamp


def internal_element(lamp: float, ff: float):
    """
    EN1993-1-5 Section 4.4(2) Equation 4.2
    :param lamp: Relative slenderness.
    :param ff: Stress ratio.
    :return: rho, Reduction factor.
    """
    rho = None
    if lamp <= 0.5 + math.sqrt(0.085 - 0.055 * ff):
        rho = 1.0
    else:
        rho = (lamp - 0.055 * (3 + ff)) / math.pow(lamp, 2)
        if rho > 1.0: rho = 1.0
    return rho


def outstand_element(lamp: float):
    """
    EN1993-1-5 Section 4.4(2) Equation 4.3
    :param lamp: Relative slenderness.
    :param ff: Stress ratio.
    :return: rho, Reduction factor.
    """
    rho = None
    if lamp <= 0.748:
        rho = 1.0
    else:
        rho = (lamp - 0.188) / math.pow(lamp, 2)
        if rho > 1.0: rho = 1.0
    return rho


def stres_ratio(bcomp: float, bwhole: float):
    return (bcomp - bwhole) / bcomp


def Table4_1_ksigma(ff: float):
    """

    :param ff: stress ratio
    :return: buckling factor
    """
    # Finding the buckling factor as per EN 1993-1-5 Table 4.1
    if ff == 1.0:
        ksigma = 4.0
    elif 1.0 > ff > 0.0:
        ksigma = 8.2 / (1.05 + ff)
    elif ff == 0.0:
        ksigma = 7.81
    elif 0.0 > ff > -1.0:
        ksigma = 7.81 - 6.29 * ff + 9.78 * ff ** 2
    elif ff == -1:
        ksigma = 23.9
    elif -1.0 > ff > -3.0:
        ksigma = 5.98 * (1.0 - ff) ** 2
    else:
        ksigma = "No Value! ksigma < -3.0"
    return ksigma


def Table4_2_ksigma(ff: float, iscompend: bool):
    """
    EN1993-1-5 Table 4.2
    :param ff: stress ratio
    :param iscompend: bool, If compression at the free end True.
    :return: buckling factor.
    """
    # Finding the buckling factor as per EN 1993-1-5 Table 4.2
    if iscompend:
        ksigma = 0.57 - 0.21 * ff + 0.07 * ff ** 2
    else:
        if ff == 1.0:
            ksigma = 0.43
        elif 1.0 > ff > 0.0:
            ksigma = 0.578 / (ff + 0.34)
        elif ff == 0.0:
            ksigma = 1.70
        elif 0.0 > ff > -1.0:
            ksigma = 1.7 - 5.0 * ff + 17.1 * ff ** 2
        elif ff == -1.0:
            ksigma = 23.8
        else:
            ksigma = "No Value! ksigma < -1.0"
    return ksigma


def Table4_1_beff(ff: float, b: float, rho: float):
    # Calculating the effective length as per EN 1993-1-5 Table 4.1
    if ff == 1.0:
        beff = rho * b
        be1 = 0.5 * beff
        be2 = 0.5 * beff
    elif 1.0 > ff >= 0.0:
        beff = rho * b
        be1 = 2 / (5 - ff) * beff
        be2 = beff - be1
    else:
        beff = rho * b / (1 - ff)
        be1 = 0.4 * beff
        be2 = 0.6 * beff
    return beff, be1, be2


def Table4_2_beff(b: float, rho: float, ff: float):
    """
    EN1993-1-5 Table 4.2
    :param ff: stress ratio
    :param b: Element width.
    :param rho: Reduction factor.
    :return: effective length, length in compression side, length in tension side.
    """
    # Calculating the effective length as per EN 1993-1-5 Table 4.2
    if 1.0 >= ff >= 0.0:
        bc = b  # Compression part of the flat width
        bt = 0.0  # Tension part of the flat width
        beff = rho * b
    else:
        bc = b / (1 - ff)  # Compression part of the flat width
        beff = rho * bc
        bt = b - bc  # Tension part of the flat width

    return beff, bc, bt


def Table4_1(b: float, rho: float, s1: float, s2: float):
    """
    EN1993-1-5 Table 4.1
    :param b: Element width.
    :param rho: Reduction factor.
    :param s1: Larger stress value on the element. It always positive (compression).
    :param s2: Equal or the smaller stress on the element. It may be negative (tension).
    :return: buckling factor, effective length.
    """
    # Finding the buckling factor as per EN 1993-1-5 Table 4.1
    ff = s2 / s1  # The stress ratio
    if ff == 1.0:
        ksigma = 4.0
    elif 1.0 > ff > 0.0:
        ksigma = 8.2 / (1.05 + ff)
    elif ff == 0.0:
        ksigma = 7.81
    elif 0.0 > ff > -1.0:
        ksigma = 7.81 - 6.29 * ff + 9.78 * ff ** 2
    elif ff == -1:
        ksigma = 23.9
    elif -1.0 > ff > -3.0:
        ksigma = 5.98 * (1.0 - ff) ** 2
    else:
        ksigma = "No Value! ksigma < -3.0"

    # Calculating the effective length as per EN 1993-1-5 Table 4.1
    if ff == 1.0:
        beff = rho * b
        be1 = 0.5 * beff
        be2 = 0.5 * beff
    elif 1.0 > ff >= 0.0:
        beff = rho * b
        be1 = 2 / (5 - ff) * beff
        be2 = beff - be1
    else:
        beff = rho * b / (1 - ff)
        be1 = 0.4 * beff
        be2 = 0.6 * beff
    return ksigma, beff, be1, be2

