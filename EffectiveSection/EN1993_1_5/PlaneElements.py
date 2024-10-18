import math


def Sec44_in(s1, s2, fy, b, t):
    """
    Ref: 1993-1-5 Section 4.4
    Plane Elements Without Longitudinal Stiffeners | Internal compression element.
    :param s2: Larger stress value on the element. It always positive (compression).
    :param s1: Equal or the smaller stress on the element. It may be negative (tension).
    :param fy: The steel yield stress
    :param ksigma: The buckling factor corresponding to the stress ratio and boundary conditions
    :param b: The flat width
    :param t: The steel thickness
    :return: rho (Reduction factor), lamP (The relative slenderness ratio)
    """
    # Finding the buckling factor as per EN 1993-1-5 Table 4.1
    ff = s2 / s1                                                    # The stress ratio
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

    eps = math.sqrt(235 / fy)                                       # Factor depending on fy
    lamP = (b / t) / (28.4 * eps * math.sqrt(ksigma))               # The relative slenderness

    # Equations in section 4.4 (2)
    if lamP <= 0.5 + math.sqrt(0.085 - 0.055 * ff):
        rho = 1.0
    else:
        rho = (lamP - 0.055 * (3 + ff)) / math.pow(lamP, 2) <= 1.0

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

    # Results
    Results = {'rho': rho,
               'ff': ff,
               'ksigma': ksigma,
               'lamP': lamP,
               'beff': beff,
               'be1': be1,
               'be2': be2}

    return Results


def Sec44_out(s1, s2, fy, b, t, iscompend):
    """
    Ref: 1993-1-5 Section 4.4
    Plane Elements Without Longitudinal Stiffeners | Outstand compression element.
    :param iscompend: bool, If compression at the free end True.
    :param s2: Larger stress value on the element. It always positive (compression).
    :param s1: Equal or the smaller stress on the element. It may be negative (tension).
    :param fy: The steel yield stress
    :param ksigma: The buckling factor corresponding to the stress ratio and boundary conditions
    :param b: The flat width
    :param t: The steel thickness
    :return: rho (Reduction factor), lamP (The relative slenderness ratio)
    """
    # Finding the buckling factor as per EN 1993-1-5 Table 4.2
    ff = s2 / s1  # The stress ratio
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

    eps = math.sqrt(235 / fy)                                   # Factor depending on fy
    lamP = (b / t) / (28.4 * eps * math.sqrt(ksigma))           # The relative slenderness

    # Equations in section 4.4 (2)
    if lamP <= 0.748:
        rho = 1.0
    else:
        rho = (lamP - 0.188) / math.pow(lamP, 2) <= 1.0

    # Calculating the effective length as per EN 1993-1-5 Table 4.2
    if 1.0 > ff >= 0.0:
        bc = b                                                  # Compression part of the flat width
        bt = 0.0                                                # Tension part of the flat width
        beff = rho * b
    else:
        bc = b / (1 - ff)                                       # Compression part of the flat width
        beff = rho * bc
        bt = b - bc                                             # Tension part of the flat width

    # Results
    Results = {'rho': rho,
               'ff': ff,
               'ksigma': ksigma,
               'lamP': lamP,
               'beff': beff,
               'bc': bc,
               'bt': bt
               }

    return Results
