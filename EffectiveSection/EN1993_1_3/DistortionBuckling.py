"""
5.5.3.2 Plane elements with edge stiffeners
"""
import math


def checkEdgeFold(b, t, doublefold: bool):
    """
    This function checks the flange width to thickness ratio is ok or not.
    :param b: Flange width.
    :param t: Thickness
    :param doublefold: Bool double edge stiffener is True
    :return: Edge stiffener is ok or not
    """
    if doublefold:
        if b / t > 90:
            isEdgeStiffener_OK = False
        else:
            isEdgeStiffener_OK = True
    else:
        if b / t > 60:
            isEdgeStiffener_OK = False
        else:
            isEdgeStiffener_OK = True
    return isEdgeStiffener_OK


def effWidthofEdgeFold(bp, be2, t, rhoc, rhod, bpc, bpd, doublefold: bool):
    """
    This function calculates the effective width of the stiffener.
    EN1993_1_3 section 5.5.3.2 Plane elements with edge stiffeners.
    :param bp: Flange width.
    :param be2: Effective flange width close to the edge stiffener.
    :param t: Thickness.
    :param rhoc: Reduction factor for the first fold of stiffener.
    :param rhod: Reduction factor for the second fold of stiffener.
    :param bpc: Length of the first fold of the stiffener.
    :param bpd: Length of the second fold of the stiffener.
    :param doublefold: Bool double edge stiffener is True.
    :return: Effective section lengths, area and buckling factor.
    """


    if doublefold:
        deff = rhod * bpd
        ceff = rhoc * bpc
        As = t * (be2 + ceff + deff)
        return [deff, ceff, As]
    else:
        ceff = rhoc * bpc
        As = t * (be2 + ceff)
        ksc = 0
        if bpc / bp <= 0.35:
            ksc = 0.5
        elif 0.35 < bpc / bp <= 0.6:
            ksc = 0.5 + 0.83 * math.pow((bpc / bp - 0.35), (2 / 3))
        else:
            print(f"bpc/bp = {bpc / bp:.3f} > 0.6, check EN1993_1_3 5.5.3.2 (5)")
        return [ceff, As, ksc]


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


def ksid(fy, sigmacrs):
    """
    This function calculates reduction factor ksid. EN1993_1_3 5.5.3.1 (7)
    :param fy: Steel yield stress
    :param sigmacrs: the elastic critical stress for the stiffeners.
    :return: The reduction factor.
    """
    lambdaD = math.sqrt(fy / sigmacrs)
    if lambdaD <= 0.65:
        ksi = 1.0
    elif 0.65 < lambdaD < 1.38:
        ksi = 1.47 - 0.723 * lambdaD
    else:
        ksi = 0.66 / lambdaD
    return ksi
