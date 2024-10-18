import math


def ksig(bpc: float, bp: float):
    """
    EN 1993-1-3 Section 5.5.3.2(5)
    :param bpc: Length of the lip.
    :param bp: Length of the flange.
    :return: Buckling factor.
    """
    ksigma = None
    if bpc / bp <= 0.35:
        ksigma = 0.5
    elif 0.35 < bpc / bp <= 0.60:
        ksigma = 0.50 + 0.83 * math.pow(bpc / bp - 0.35, (2.0 / 3.0))
    else:
        print(f"bpc/bp = {bpc / bp:.3f} > 0.6, check EN1993_1_3 5.5.3.2 (5)")
    return ksigma


def calc_b1(b, be2, t, ceff):
    return b - (be2 * t * be2 / 2) / ((be2 + ceff) * t)


def springStiffnessK(E, t, v, b1, hw, b2, bending: bool):
    """
    This function calculates the spring stiffness as per EN1993_1_3 5.5.3.1(5)
    :param E: Elastic modulus
    :param t: Thickness
    :param v: Poisson's ratio
    :param b1: Distance
    :param hw: Distance
    :param b2: Distance
    :param bending: Section is subjected to bending is True
    :return: K spring stiffness
    """
    if bending:
        kf = 0
    else:
        kf = 1
    K = (E * t ** 3) / (4 * (1 - v ** 2)) * (1 / (b1 ** 2 * hw + b1 ** 3 + 0.5 * b1 * b2 * hw * kf))
    return K


def Is(be2, t, ceff):
    Is = be2 * t ** 3 / 12 + ceff ** 3 * t / 12 + be2 * t * (ceff ** 2 / (2 * (be2 + ceff))) ** 2 + ceff * t * (
            ceff / 2 - ceff ** 2 / (2 * (be2 + ceff))) ** 2
    return Is


def calc_scrs(K, Is, E, As):
    scrs = 2 * math.sqrt(K * E * Is) / As
    return scrs


def thk_reduction(fy, scrs):
    lamd = math.sqrt(fy / scrs)
    if lamd <= 0.65:
        xd = 1.0
    elif 0.65 < lamd <= 1.38:
        xd = 1.47 - 0.723 * lamd
    else:
        xd = 0.66 / lamd
    return xd

