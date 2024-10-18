import math
import matplotlib.pyplot as plt
import numpy as np
import EffectiveSection.EN1993_1_5.Sec4 as Sec4
import EffectiveSection.EN1993_1_3.Sec5_5_3 as Sec553
import EffectiveSection.Modes.IntCalcProp as intprop


class GrossProp:
    def __init__(self, a: float, b: float, c: float, t: float, ro: float, f: float):
        self.elements = None
        self.nodes = None
        self.y_inches = None
        self.x_inches = None
        self.descp = None
        self.y = None
        self.x = None
        self.tcore = None
        self.cc = None
        self.bb = None
        self.aa = None
        self.r = None
        self.c = None
        self.b = None
        self.a = None
        self.R = ro
        self.t = t
        self.C = c
        self.B = b
        self.A = a
        self.xo = None
        self.It = None
        self.Cw = None
        self.ysc = None
        self.xsc = None
        self.Iw = None
        self.Ixy = None
        self.Wy = None
        self.Iy = None
        self.Wx = None
        self.Ix = None
        self.zgy = None
        self.zgx = None
        self.Ar = None
        self.E = 210000.0  # MPa
        self.v = 0.3  # Poisson's ratio

        # Section parts
        #         6       7
        #       ┌────   ────┐    ↑
        #       │           │ 8  │
        #       │ 5              hc
        #                        │
        #       │                ↓
        #     --│-------------------
        #       │
        #       │ 4
        #       │           │ 1
        #       └────   ────┘
        #          3      2

        secDivider = '============================================'
        Rep = f'{secDivider}\nUSER INPUT\n{secDivider}\n'
        Rep += f'==== Sections ====\n'
        Rep += f'A = {self.A:.2f} mm\n'
        Rep += f'B = {self.B:.2f} mm\n'
        Rep += f'C = {self.C:.2f} mm\n'
        Rep += f't = {self.t:.2f} mm\n'
        Rep += f'R = {self.R:.2f} mm\n'
        print(Rep)
        self.lippedCSection()
        self.grossProp(self.x, self.y, self.t, self.r)

        # Variables for axial compression
        self.scomed = f
        self.Axial_elementData2 = None
        self.Axial_Aeff = None
        self.Axial_ygct = None
        self.Axial_ygc = None
        self.Axial_xgc = None
        self.Axial_dxgc = None
        self.calcs_AxialCompression()
        print(self.Axial_Aeff)
        # Variables for bending about strong axis
        self.BendStrong_elementData2 = None
        self.BendStrong_ygct = None  # Top flange extreme fiber to neutral axis
        self.BendStrong_ygc = None  # Bottom flange extreme fiber to neutral axis
        self.BendStrong_Ixeff = None  # Moment of inertia
        self.BendStrong_Wxeff = None  # Section modulus
        self.calcs_BendingStrong()
        # Variables for bending about strong axis
        self.BendWeakLip_xgct = None
        self.BendWeakLip_xgc = None
        self.BendWeakLip_elementData2 = None
        self.BendWeakLip_ygct = None  # Top flange extreme fiber to neutral axis
        self.BendWeakLip_ygc = None  # Bottom flange extreme fiber to neutral axis
        self.BendWeakLip_Iyeff = None  # Moment of inertia
        self.BendWeakLip_Wyeff = None  # Section modulus
        self.calcs_BendingWeakLip()
        # Variables for bending about strong axis
        self.BendWeakWeb_xgct = None
        self.BendWeakWeb_xgc = None
        self.BendWeakWeb_elementData2 = None
        self.BendWeakWeb_ygct = None  # Top flange extreme fiber to neutral axis
        self.BendWeakWeb_ygc = None  # Bottom flange extreme fiber to neutral axis
        self.BendWeakWeb_Iyeff = None  # Moment of inertia
        self.BendWeakWeb_Wyeff = None  # Section modulus
        self.calcs_BendingWeakWeb()
        self.plot_effC_section()

    def lippedCSection(self):
        r = self.R + self.t / 2.0
        # Centerline dimensions
        aa = self.A - self.t
        bb = self.B - self.t
        cc = self.C - self.t / 2.0
        tcore = self.t - 0.04
        # Flat portions
        a = aa - 2 * r
        b = bb - 2 * r
        c = cc - r

        self.a = a
        self.b = b
        self.c = c
        self.r = r
        self.aa = aa
        self.bb = bb
        self.cc = cc
        self.tcore = tcore

        def tranform(alfa, origin):
            # Transformation matrix
            p = origin
            c, s = np.cos(math.radians(alfa)), np.sin(math.radians(alfa))
            j = np.array([[c, s, 0],
                          [-s, c, 0],
                          [0, 0, 1]])
            RotatedCorners = np.matmul(j, p)
            return RotatedCorners.T

        # Bottom right
        origin1 = np.array([bb - r, r, 0.0])
        radius = np.array([0, r, 0.0])
        start_ang = 90
        p11 = tranform(start_ang + 10, radius)
        p12 = tranform(start_ang + 20, radius)
        p13 = tranform(start_ang + 30, radius)
        p14 = tranform(start_ang + 40, radius)
        p15 = tranform(start_ang + 50, radius)
        p16 = tranform(start_ang + 60, radius)
        p17 = tranform(start_ang + 70, radius)
        p18 = tranform(start_ang + 80, radius)
        # Bottom left
        origin2 = np.array([r, r, 0.0])
        start_ang = 180
        p21 = tranform(start_ang + 10, radius)
        p22 = tranform(start_ang + 20, radius)
        p23 = tranform(start_ang + 30, radius)
        p24 = tranform(start_ang + 40, radius)
        p25 = tranform(start_ang + 50, radius)
        p26 = tranform(start_ang + 60, radius)
        p27 = tranform(start_ang + 70, radius)
        p28 = tranform(start_ang + 80, radius)
        # Top left
        origin3 = np.array([r, aa - r, 0.0])
        start_ang = 270
        p31 = tranform(start_ang + 10, radius)
        p32 = tranform(start_ang + 20, radius)
        p33 = tranform(start_ang + 30, radius)
        p34 = tranform(start_ang + 40, radius)
        p35 = tranform(start_ang + 50, radius)
        p36 = tranform(start_ang + 60, radius)
        p37 = tranform(start_ang + 70, radius)
        p38 = tranform(start_ang + 80, radius)
        # Top right
        origin4 = np.array([bb - r, aa - r, 0.0])
        start_ang = 0
        p41 = tranform(start_ang + 10, radius)
        p42 = tranform(start_ang + 20, radius)
        p43 = tranform(start_ang + 30, radius)
        p44 = tranform(start_ang + 40, radius)
        p45 = tranform(start_ang + 50, radius)
        p46 = tranform(start_ang + 60, radius)
        p47 = tranform(start_ang + 70, radius)
        p48 = tranform(start_ang + 80, radius)

        self.nodes = np.array([[0, bb, cc, 1, 1, 1, 1, 0],
                               [1, bb, r, 1, 1, 1, 1, 0],
                               [2, origin1[0] + p11[0], origin1[1] + p11[1], 1, 1, 1, 1, 0],
                               [3, origin1[0] + p12[0], origin1[1] + p12[1], 1, 1, 1, 1, 0],
                               [4, origin1[0] + p13[0], origin1[1] + p13[1], 1, 1, 1, 1, 0],
                               [5, origin1[0] + p14[0], origin1[1] + p14[1], 1, 1, 1, 1, 0],
                               [6, origin1[0] + p15[0], origin1[1] + p15[1], 1, 1, 1, 1, 0],
                               [7, origin1[0] + p16[0], origin1[1] + p16[1], 1, 1, 1, 1, 0],
                               [8, origin1[0] + p17[0], origin1[1] + p17[1], 1, 1, 1, 1, 0],
                               [9, origin1[0] + p18[0], origin1[1] + p18[1], 1, 1, 1, 1, 0],
                               [10, bb - r, 0, 1, 1, 1, 1, 0],
                               [11, r + b / 2.0, 0, 1, 1, 1, 1, 0],
                               [12, r, 0, 1, 1, 1, 1, 0],
                               [13, origin2[0] + p21[0], origin2[1] + p21[1], 1, 1, 1, 1, 0],
                               [14, origin2[0] + p22[0], origin2[1] + p22[1], 1, 1, 1, 1, 0],
                               [15, origin2[0] + p23[0], origin2[1] + p23[1], 1, 1, 1, 1, 0],
                               [16, origin2[0] + p24[0], origin2[1] + p24[1], 1, 1, 1, 1, 0],
                               [17, origin2[0] + p25[0], origin2[1] + p25[1], 1, 1, 1, 1, 0],
                               [18, origin2[0] + p26[0], origin2[1] + p26[1], 1, 1, 1, 1, 0],
                               [19, origin2[0] + p27[0], origin2[1] + p27[1], 1, 1, 1, 1, 0],
                               [20, origin2[0] + p28[0], origin2[1] + p28[1], 1, 1, 1, 1, 0],
                               [21, 0, r, 1, 1, 1, 1, 0],
                               [22, 0, r + a * (1.0 / 4.0), 1, 1, 1, 1, 0],
                               [23, 0, r + a * (2.0 / 4.0), 1, 1, 1, 1, 0],
                               [24, 0, r + a * (3.0 / 4.0), 1, 1, 1, 1, 0],
                               [25, 0, r + a, 1, 1, 1, 1, 0],
                               [26, origin3[0] + p31[0], origin3[1] + p31[1], 1, 1, 1, 1, 0],
                               [27, origin3[0] + p32[0], origin3[1] + p32[1], 1, 1, 1, 1, 0],
                               [28, origin3[0] + p33[0], origin3[1] + p33[1], 1, 1, 1, 1, 0],
                               [29, origin3[0] + p34[0], origin3[1] + p34[1], 1, 1, 1, 1, 0],
                               [30, origin3[0] + p35[0], origin3[1] + p35[1], 1, 1, 1, 1, 0],
                               [31, origin3[0] + p36[0], origin3[1] + p36[1], 1, 1, 1, 1, 0],
                               [32, origin3[0] + p37[0], origin3[1] + p37[1], 1, 1, 1, 1, 0],
                               [33, origin3[0] + p38[0], origin3[1] + p38[1], 1, 1, 1, 1, 0],
                               [34, r, aa, 1, 1, 1, 1, 0],
                               [35, r + b / 2.0, aa, 1, 1, 1, 1, 0],
                               [36, bb - r, aa, 1, 1, 1, 1, 0],
                               [37, origin4[0] + p41[0], origin4[1] + p41[1], 1, 1, 1, 1, 0],
                               [38, origin4[0] + p42[0], origin4[1] + p42[1], 1, 1, 1, 1, 0],
                               [39, origin4[0] + p43[0], origin4[1] + p43[1], 1, 1, 1, 1, 0],
                               [40, origin4[0] + p44[0], origin4[1] + p44[1], 1, 1, 1, 1, 0],
                               [41, origin4[0] + p45[0], origin4[1] + p45[1], 1, 1, 1, 1, 0],
                               [42, origin4[0] + p46[0], origin4[1] + p46[1], 1, 1, 1, 1, 0],
                               [43, origin4[0] + p47[0], origin4[1] + p47[1], 1, 1, 1, 1, 0],
                               [44, origin4[0] + p48[0], origin4[1] + p48[1], 1, 1, 1, 1, 0],
                               [45, bb, aa - r, 1, 1, 1, 1, 0],
                               [46, bb, aa - cc, 1, 1, 1, 1, 0]])

        self.elements = np.array([[0, 0, 1, self.t, 0],
                                  [1, 1, 2, self.t, 0],
                                  [2, 2, 3, self.t, 0],
                                  [3, 3, 4, self.t, 0],
                                  [4, 4, 5, self.t, 0],
                                  [5, 5, 6, self.t, 0],
                                  [6, 6, 7, self.t, 0],
                                  [7, 7, 8, self.t, 0],
                                  [8, 8, 9, self.t, 0],
                                  [9, 9, 10, self.t, 0]
                                  ])
        self.x = self.nodes[:, 1]
        self.y = self.nodes[:, 2]
        self.descp = f'Section : A: {self.A:.2f}, B: {self.B:.2f}, C: {self.C:.2f}, t: {self.t:.2f}'
        self.x_inches = self.x * 25.4
        self.y_inches = self.y * 25.4
        # Data for plotting
        self.plot_C_section()

    def grossProp(self, x, y, t, r):
        # Area of cross-section
        da = np.zeros([len(x)])
        ba = np.zeros([len(x)])
        for i in range(1, len(da)):
            da[i] = math.sqrt(math.pow(x[i - 1] - x[i], 2) + math.pow(y[i - 1] - y[i], 2)) * t
            ba[i] = math.sqrt(math.pow(x[i - 1] - x[i], 2) + math.pow(y[i - 1] - y[i], 2))
        Ar = np.sum(da)
        Lt = np.sum(ba)
        # Total rj.tetaj/90
        Trj = 4 * 4 * r * (1.0 / 4.0)
        delta = 0.43 * Trj / Lt
        # First moment of area and coordinate for gravity centre
        sx0 = np.zeros([len(x)])
        sy0 = np.zeros([len(x)])
        for i in range(1, len(sx0)):
            sx0[i] = (y[i] + y[i - 1]) * da[i] / 2
        zgy = np.sum(sx0) / Ar
        for i in range(1, len(sy0)):
            sy0[i] = (x[i] + x[i - 1]) * da[i] / 2
        zgx = np.sum(sy0) / Ar

        # Second moment of area
        Ix0 = np.zeros([len(x)])
        Iy0 = np.zeros([len(x)])
        for i in range(1, len(Ix0)):
            Ix0[i] = (math.pow(y[i], 2) + math.pow(y[i - 1], 2) + y[i] * y[i - 1]) * da[i] / 3
        for i in range(1, len(Iy0)):
            Iy0[i] = (math.pow(x[i], 2) + math.pow(x[i - 1], 2) + x[i] * x[i - 1]) * da[i] / 3
        Ix = np.sum(Ix0) - Ar * math.pow(zgy, 2)
        Iy = np.sum(Iy0) - Ar * math.pow(zgx, 2)

        # Product moment of area
        Ixy0 = np.zeros([len(x)])
        for i in range(1, len(Ixy0)):
            Ixy0[i] = (2 * x[i - 1] * y[i - 1] + 2 * x[i] * y[i] + x[i - 1] * y[i] + x[i] * y[i - 1]) * da[i] / 6
        Ixy = np.sum(Ixy0) - (np.sum(sx0) * np.sum(sy0)) / Ar

        # Principle axis
        alfa = 0.5 * math.atan(2 * Ixy / (Iy - Ix))
        Iksi = 0.5 * (Ix + Iy + math.sqrt(math.pow(Iy - Ix, 2) + 4 * math.pow(Ixy, 2)))
        Ieta = 0.5 * (Ix + Iy - math.sqrt(math.pow(Iy - Ix, 2) + 4 * math.pow(Ixy, 2)))

        # Sectoral coordinates
        w = np.zeros([len(x)])
        w0 = np.zeros([len(x)])
        Iw = np.zeros([len(x)])
        w0[0] = 0
        for i in range(1, len(w0)):
            w0[i] = x[i - 1] * y[i] - x[i] * y[i - 1]
            w[i] = w[i - 1] + w0[i]
            Iw[i] = (w[i - 1] + w[i]) * da[i] / 2
        wmean = np.sum(Iw) / Ar

        # Sectorial constants
        Ixw0 = np.zeros([len(x)])
        Iyw0 = np.zeros([len(x)])
        Iww0 = np.zeros([len(x)])
        for i in range(1, len(Ixw0)):
            Ixw0[i] = (2 * x[i - 1] * w[i - 1] + 2 * x[i] * w[i] + x[i - 1] * w[i] + x[i] * w[i - 1]) * da[i] / 6
            Iyw0[i] = (2 * y[i - 1] * w[i - 1] + 2 * y[i] * w[i] + y[i - 1] * w[i] + y[i] * w[i - 1]) * da[i] / 6
            Iww0[i] = (math.pow(w[i], 2) + math.pow(w[i - 1], 2) + w[i] * w[i - 1]) * da[i] / 3
        Ixw = np.sum(Ixw0) - np.sum(sy0) * np.sum(Iw) / Ar
        Iyw = np.sum(Iyw0) - np.sum(sx0) * np.sum(Iw) / Ar
        Iww = np.sum(Iww0) - math.pow(np.sum(Iw), 2) / Ar

        # Shear centre
        xsc = (Iyw * Iy - Ixw * Ixy) / (Ix * Iy - math.pow(Ixy, 2))
        ysc = (-Ixw * Ix + Iyw * Ixy) / (Ix * Iy - math.pow(Ixy, 2))

        # Warping constant
        Cw = Iww + ysc * Ixw - xsc * Iyw

        # Torsion constant
        It0 = np.zeros([len(x)])
        for i in range(1, len(It0)):
            It0[i] = da[i] * math.pow(t, 2) / 3
        It = np.sum(It0)

        # Distance between centroid and shear centre
        xo = abs(xsc) + zgx
        # Distances from the boundaries
        zgb = zgy
        zgt = max(y) - zgb
        zgl = zgx
        zgr = max(x) - zgl
        # Calculated data
        propData = np.zeros([14])
        propData[0] = Ar * (1 - delta)
        propData[1] = max(zgb, zgt)
        propData[2] = max(zgl, zgr)
        propData[3] = Ix * (1 - 2 * delta)
        propData[4] = Ix * (1 - 2 * delta) / max(zgb, zgt)
        propData[5] = Iy * (1 - 2 * delta)
        propData[6] = Iy * (1 - 2 * delta) / max(zgl, zgr)
        propData[7] = Ixy
        propData[8] = np.sum(Iw)
        propData[9] = xsc
        propData[10] = ysc
        propData[11] = Cw * (1 - 4 * delta)
        propData[12] = It * (1 - 2 * delta)
        propData[13] = xo
        Wx = Ix * (1 - 2 * delta) / max(zgb, zgt)
        Wy = Iy * (1 - 2 * delta) / max(zgl, zgr)
        # Define the attributes
        self.Ar = Ar
        self.zgx = zgx
        self.zgy = zgy
        self.Ix = Ix
        self.Wx = Wx
        self.Iy = Iy
        self.Wy = Wy
        self.Ixy = Ixy
        self.Iw = Iw
        self.xsc = xsc
        self.ysc = ysc
        self.Cw = Cw
        self.It = It
        self.xo = xo

        # Data dictionary
        self.prop = {
            "Ag": [Ar, ' mm2, Area of cross-section'],
            "zgx": [zgx, ' mm, Coordinate for gravity centre'],
            "zgy": [zgy, ' mm, Coordinate for gravity centre'],
            "Ix": [Ix, ' mm4, Second moment of area about strong axis'],
            "Wx": [Ix * (1 - 2 * delta) / max(zgb, zgt), ' mm3, Section modulus about strong axis'],
            "Iy": [Iy, ' mm4, Second moment of area about weak axis'],
            "Wy": [Iy * (1 - 2 * delta) / max(zgl, zgr), ' mm3, Section modulus about weak axis'],
            "Ixy": [Ixy, ' mm4, Product moment of area'],
            "Iw": [np.sum(Iw), ' Sectorial constant'],
            "xsc": [xsc, ' mm, Shear center on y axis'],
            "ysc": [ysc, ' mm, Shear center on x axis'],
            "Cw": [Cw, ' mm6, Warping constant'],
            "It": [It, ' mm4, Torsional constant'],
            "xo": [xo, ' mm, Distance between centroid and shear centre']
        }
        # print(self.prop)
        secDivider = '============================================'
        print(f'{secDivider}\nGROSS SECTION PROPERTIES (in mm)\n{secDivider}')
        for key, value in self.prop.items():
            print(f'{key}: {value[0]:.4f}{value[1]}')

        return self.prop, propData

    def calcs_AxialCompression(self):
        # Design Stress
        scomed = self.scomed
        # ==============================================================================================================
        # Effective width of the top flange
        # ==============================================================================================================
        top_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        top_flg_ksigma = Sec4.Table4_1_ksigma(top_flg_stress_ratio)
        top_flg_lamp = Sec4.lamp(self.bb, self.tcore, top_flg_ksigma, scomed, True)
        top_flg_rho = Sec4.internal_element(top_flg_lamp, top_flg_stress_ratio)
        top_flg_beff = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[0]
        top_flg_be1 = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[1]
        top_flg_be2 = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[2]
        # ==============================================================================================================
        # Effective width of the bot flange
        # ==============================================================================================================
        bot_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        bot_flg_ksigma = Sec4.Table4_1_ksigma(bot_flg_stress_ratio)
        bot_flg_lamp = Sec4.lamp(self.bb, self.tcore, bot_flg_ksigma, scomed, True)
        bot_flg_rho = Sec4.internal_element(bot_flg_lamp, bot_flg_stress_ratio)
        bot_flg_beff = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[0]
        bot_flg_be1 = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[1]
        bot_flg_be2 = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[2]
        # ==============================================================================================================
        # Effective width of the top edge fold
        # ==============================================================================================================
        top_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        top_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        top_lip_lamp = Sec4.lamp(self.cc, self.tcore, top_lip_ksigma, scomed, True)
        top_lip_rho = Sec4.outstand_element(top_lip_lamp)
        top_lip_beff = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[0]
        top_lip_bc = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[1]
        top_lip_bt = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        top_As = self.tcore * (top_flg_be2 + top_lip_beff)
        # Spring stiffness of the edge stiffener
        top_b1 = Sec553.calc_b1(self.bb, top_flg_be2, self.tcore, top_lip_beff)
        top_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, top_b1, self.aa,
                                        top_b1, False)
        top_Is = Sec553.Is(top_flg_be2, self.tcore, top_lip_beff)
        top_scrs = Sec553.calc_scrs(top_K, top_Is, self.E, top_As)
        # Thickness reduction factor
        top_xd = Sec553.thk_reduction(scomed, top_scrs)
        top_t_red = top_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the bottom edge fold
        # ==============================================================================================================
        bot_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        bot_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        bot_lip_lamp = Sec4.lamp(self.cc, self.tcore, bot_lip_ksigma, scomed, True)
        bot_lip_rho = Sec4.outstand_element(bot_lip_lamp)
        bot_lip_beff = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[0]
        bot_lip_bc = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[1]
        bot_lip_bt = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        bot_As = self.tcore * (bot_flg_be2 + bot_lip_beff)
        # Spring stiffness of the edge stiffener
        bot_b1 = Sec553.calc_b1(self.bb, bot_flg_be2, self.tcore, bot_lip_beff)
        bot_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, bot_b1, self.aa,
                                        bot_b1, False)
        bot_Is = Sec553.Is(bot_flg_be2, self.tcore, bot_lip_beff)
        bot_scrs = Sec553.calc_scrs(bot_K, bot_Is, self.E, bot_As)
        # Thickness reduction factor
        bot_xd = Sec553.thk_reduction(scomed, bot_scrs)
        bot_t_red = bot_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the web
        # ==============================================================================================================
        elementData1 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])
        # Distance from tension (bottom) flange
        ht = intprop.calcProps(elementData1)[1]
        # Distance from compression (top) flange
        hc = self.aa - ht
        # Stress ratio
        ff = 1.0  # uniform compression
        web_ksigma = Sec4.Table4_1_ksigma(ff)
        web_lamp = Sec4.lamp(self.aa, self.tcore, web_ksigma, scomed, True)
        web_rho = Sec4.internal_element(web_lamp, ff)
        web_beff = Sec4.Table4_1_beff(ff, self.aa, web_rho)[0]
        web_be1 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[1]
        web_be2 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[2]
        h1 = web_be1
        h2 = web_be2
        # ==============================================================================================================
        # Effective section properties
        # ==============================================================================================================
        # Create a matrix contains the element data from bottom lip to top lip
        # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
        self.Axial_elementData2 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [4, 0.0, 0.0, 0.0, h2, self.tcore],
             [5, 0.0, self.aa - h1, 0.0, self.aa, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])

        # Results
        self.Axial_ygc = intprop.calcProps(self.Axial_elementData2)[1]
        self.Axial_ygct = self.aa - intprop.calcProps(self.Axial_elementData2)[1]
        self.Axial_xgc = intprop.calcProps(self.Axial_elementData2)[3]
        self.Axial_Aeff = intprop.calcProps(self.Axial_elementData2)[0]
        self.Axial_dxgc = self.Axial_xgc - self.zgx  # if it is + compression on web.

    def calcs_BendingStrong(self):
        # Design Stress
        scomed = self.scomed
        # ==============================================================================================================
        # Effective width of the top flange
        # ==============================================================================================================
        top_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        top_flg_ksigma = Sec4.Table4_1_ksigma(top_flg_stress_ratio)
        top_flg_lamp = Sec4.lamp(self.bb, self.tcore, top_flg_ksigma, scomed, True)
        top_flg_rho = Sec4.internal_element(top_flg_lamp, top_flg_stress_ratio)
        top_flg_beff = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[0]
        top_flg_be1 = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[1]
        top_flg_be2 = Sec4.Table4_1_beff(top_flg_stress_ratio, self.bb, top_flg_rho)[2]
        # ==============================================================================================================
        # Effective width of the bot flange
        # ==============================================================================================================
        bot_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        bot_flg_ksigma = Sec4.Table4_1_ksigma(bot_flg_stress_ratio)
        bot_flg_lamp = Sec4.lamp(self.bb, self.tcore, bot_flg_ksigma, scomed, False)
        bot_flg_rho = Sec4.internal_element(bot_flg_lamp, bot_flg_stress_ratio)
        bot_flg_beff = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[0]
        bot_flg_be1 = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[1]
        bot_flg_be2 = Sec4.Table4_1_beff(bot_flg_stress_ratio, self.bb, bot_flg_rho)[2]
        # ==============================================================================================================
        # Effective width of the top edge fold
        # ==============================================================================================================
        top_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        top_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        top_lip_lamp = Sec4.lamp(self.cc, self.tcore, top_lip_ksigma, scomed, True)
        top_lip_rho = Sec4.outstand_element(top_lip_lamp)
        top_lip_beff = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[0]
        top_lip_bc = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[1]
        top_lip_bt = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        top_As = self.tcore * (top_flg_be2 + top_lip_beff)
        # Spring stiffness of the edge stiffener
        top_b1 = Sec553.calc_b1(self.bb, top_flg_be2, self.tcore, top_lip_beff)
        top_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, top_b1, self.aa,
                                        top_b1, True)
        top_Is = Sec553.Is(top_flg_be2, self.tcore, top_lip_beff)
        top_scrs = Sec553.calc_scrs(top_K, top_Is, self.E, top_As)
        # Thickness reduction factor
        top_xd = Sec553.thk_reduction(scomed, top_scrs)
        top_t_red = top_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the bottom edge fold
        # ==============================================================================================================
        bot_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        bot_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        bot_lip_lamp = Sec4.lamp(self.cc, self.tcore, bot_lip_ksigma, scomed, False)
        bot_lip_rho = Sec4.outstand_element(bot_lip_lamp)
        bot_lip_beff = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[0]
        bot_lip_bc = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[1]
        bot_lip_bt = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        bot_As = self.tcore * (bot_flg_be2 + bot_lip_beff)
        # Spring stiffness of the edge stiffener
        bot_b1 = Sec553.calc_b1(self.bb, bot_flg_be2, self.tcore, bot_lip_beff)
        bot_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, bot_b1, self.aa,
                                        bot_b1, True)
        bot_Is = Sec553.Is(bot_flg_be2, self.tcore, bot_lip_beff)
        bot_scrs = Sec553.calc_scrs(bot_K, bot_Is, self.E, bot_As)
        # Thickness reduction factor
        bot_xd = 1.0  # Tension on bottom flange
        bot_t_red = bot_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the web
        # ==============================================================================================================
        elementData1 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])
        # Distance from tension (bottom) flange
        ht = intprop.calcProps(elementData1)[1]
        # Distance from compression (top) flange
        hc = self.aa - ht
        # Stress ratio
        ff = (hc - self.aa) / hc
        web_ksigma = Sec4.Table4_1_ksigma(ff)
        web_lamp = Sec4.lamp(self.aa, self.tcore, web_ksigma, scomed, True)
        web_rho = Sec4.internal_element(web_lamp, ff)
        web_beff = Sec4.Table4_1_beff(ff, self.aa, web_rho)[0]
        web_be1 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[1]
        web_be2 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[2]
        h1 = web_be1
        h2 = self.aa - (hc - web_be2)
        # ==============================================================================================================
        # Effective section properties
        # ==============================================================================================================
        # Create a matrix contains the element data from bottom lip to top lip
        # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
        self.BendStrong_elementData2 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [4, 0.0, 0.0, 0.0, h2, self.tcore],
             [5, 0.0, self.aa - h1, 0.0, self.aa, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])

        # Results
        self.BendStrong_Ixeff = intprop.calcProps(self.BendStrong_elementData2)[2]
        self.BendStrong_ygc = intprop.calcProps(self.BendStrong_elementData2)[1]
        self.BendStrong_ygct = self.aa - intprop.calcProps(self.BendStrong_elementData2)[1]
        self.BendStrong_Wxeff = intprop.calcProps(self.BendStrong_elementData2)[2] / (
                self.aa - intprop.calcProps(self.BendStrong_elementData2)[1])

    def calcs_BendingWeakLip(self):
        # Design Stress
        scomed = self.scomed
        # ==============================================================================================================
        # Effective width of the top flange
        # ==============================================================================================================
        top_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        top_flg_ksigma = Sec4.Table4_2_ksigma(top_flg_stress_ratio, True)
        top_flg_lamp = Sec4.lamp(self.bb, self.tcore, top_flg_ksigma, scomed, True)
        top_flg_rho = Sec4.internal_element(top_flg_lamp, top_flg_stress_ratio)
        top_flg_beff = Sec4.Table4_2_beff(self.bb, top_flg_rho, top_flg_stress_ratio)[0]
        top_flg_be1 = self.zgx
        top_flg_be2 = top_flg_beff
        # ==============================================================================================================
        # Effective width of the bot flange
        # ==============================================================================================================
        bot_flg_stress_ratio = Sec4.stres_ratio(self.bb, 0.0)  # bwhole is 0 for uniform stress.
        bot_flg_ksigma = Sec4.Table4_2_ksigma(bot_flg_stress_ratio, True)
        bot_flg_lamp = Sec4.lamp(self.bb, self.tcore, bot_flg_ksigma, scomed, True)
        bot_flg_rho = Sec4.internal_element(bot_flg_lamp, bot_flg_stress_ratio)
        bot_flg_beff = Sec4.Table4_2_beff(self.bb, bot_flg_rho, bot_flg_stress_ratio)[0]
        bot_flg_be1 = self.zgx
        bot_flg_be2 = bot_flg_beff
        # ==============================================================================================================
        # Effective width of the top edge fold
        # ==============================================================================================================
        top_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        top_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        top_lip_lamp = Sec4.lamp(self.cc, self.tcore, top_lip_ksigma, scomed, True)
        top_lip_rho = Sec4.outstand_element(top_lip_lamp)
        top_lip_beff = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[0]
        top_lip_bc = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[1]
        top_lip_bt = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        top_As = self.tcore * (top_flg_be2 + top_lip_beff)
        # Spring stiffness of the edge stiffener
        top_b1 = Sec553.calc_b1(self.bb, top_flg_be2, self.tcore, top_lip_beff)
        top_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, top_b1, self.aa,
                                        top_b1, False)
        top_Is = Sec553.Is(top_flg_be2, self.tcore, top_lip_beff)
        top_scrs = Sec553.calc_scrs(top_K, top_Is, self.E, top_As)
        # Thickness reduction factor
        top_xd = Sec553.thk_reduction(scomed, top_scrs)
        top_t_red = top_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the bottom edge fold
        # ==============================================================================================================
        bot_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        bot_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        bot_lip_lamp = Sec4.lamp(self.cc, self.tcore, bot_lip_ksigma, scomed, True)
        bot_lip_rho = Sec4.outstand_element(bot_lip_lamp)
        bot_lip_beff = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[0]
        bot_lip_bc = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[1]
        bot_lip_bt = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        bot_As = self.tcore * (bot_flg_be2 + bot_lip_beff)
        # Spring stiffness of the edge stiffener
        bot_b1 = Sec553.calc_b1(self.bb, bot_flg_be2, self.tcore, bot_lip_beff)
        bot_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, bot_b1, self.aa,
                                        bot_b1, False)
        bot_Is = Sec553.Is(bot_flg_be2, self.tcore, bot_lip_beff)
        bot_scrs = Sec553.calc_scrs(bot_K, bot_Is, self.E, bot_As)
        # Thickness reduction factor
        bot_xd = Sec553.thk_reduction(scomed, bot_scrs)
        bot_t_red = bot_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the web
        # ==============================================================================================================
        # Web is fully effective under tension
        ff = 1.0
        web_ksigma = Sec4.Table4_1_ksigma(ff)
        web_lamp = Sec4.lamp(self.aa, self.tcore, web_ksigma, scomed, False)
        web_rho = Sec4.internal_element(web_lamp, ff)
        web_beff = Sec4.Table4_1_beff(ff, self.aa, web_rho)[0]
        web_be1 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[1]
        web_be2 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[2]
        h1 = web_be1
        h2 = web_be1
        # ==============================================================================================================
        # Effective section properties
        # ==============================================================================================================
        # Create a matrix contains the element data from bottom lip to top lip
        # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
        self.BendWeakLip_elementData2 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [4, 0.0, 0.0, 0.0, h2, self.tcore],
             [5, 0.0, self.aa - h1, 0.0, self.aa, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])

        # Results
        self.BendWeakLip_Iyeff = intprop.calcProps(self.BendWeakLip_elementData2)[4]
        self.BendWeakLip_xgc = intprop.calcProps(self.BendWeakLip_elementData2)[3]
        self.BendWeakLip_xgct = self.bb - intprop.calcProps(self.BendWeakLip_elementData2)[3]
        self.BendWeakLip_Wyeff = intprop.calcProps(self.BendWeakLip_elementData2)[4] / (
            max(self.BendWeakLip_xgc, self.BendWeakLip_xgct))

    def calcs_BendingWeakWeb(self):
        # Design Stress
        scomed = self.scomed
        # ==============================================================================================================
        # Effective width of the top flange
        # ==============================================================================================================
        # Only tension part is taken into account for the flange (bb-zgx)
        top_flg_be1 = 0.0
        top_flg_be2 = (self.bb - self.zgx)
        # ==============================================================================================================
        # Effective width of the bot flange
        # ==============================================================================================================
        # Only tension part is taken into account for the flange (bb-zgx)
        bot_flg_be1 = 0.0
        bot_flg_be2 = (self.bb - self.zgx)
        # ==============================================================================================================
        # Effective width of the top edge fold
        # ==============================================================================================================
        # Lip is fully effective under tension.
        top_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        top_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        top_lip_lamp = Sec4.lamp(self.cc, self.tcore, top_lip_ksigma, scomed, False)
        top_lip_rho = Sec4.outstand_element(top_lip_lamp)
        top_lip_beff = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[0]
        top_lip_bc = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[1]
        top_lip_bt = Sec4.Table4_2_beff(self.cc, top_lip_rho, top_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        top_As = self.tcore * (top_flg_be2 + top_lip_beff)
        # Spring stiffness of the edge stiffener
        top_b1 = Sec553.calc_b1(self.bb, top_flg_be2, self.tcore, top_lip_beff)
        top_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, top_b1, self.aa,
                                        top_b1, True)
        top_Is = Sec553.Is(top_flg_be2, self.tcore, top_lip_beff)
        top_scrs = Sec553.calc_scrs(top_K, top_Is, self.E, top_As)
        # Thickness reduction factor
        top_xd = Sec553.thk_reduction(scomed, top_scrs)
        top_t_red = top_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the bottom edge fold
        # ==============================================================================================================
        # Lip is fully effective under tension.
        bot_lip_stres_ratio = Sec4.stres_ratio(self.cc, 0.0)  # bwhole is 0 for uniform stress.
        bot_lip_ksigma = Sec553.ksig(self.cc, self.bb)
        bot_lip_lamp = Sec4.lamp(self.cc, self.tcore, bot_lip_ksigma, scomed, False)
        bot_lip_rho = Sec4.outstand_element(bot_lip_lamp)
        bot_lip_beff = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[0]
        bot_lip_bc = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[1]
        bot_lip_bt = Sec4.Table4_2_beff(self.cc, bot_lip_rho, bot_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        bot_As = self.tcore * (bot_flg_be2 + bot_lip_beff)
        # Spring stiffness of the edge stiffener
        bot_b1 = Sec553.calc_b1(self.bb, bot_flg_be2, self.tcore, bot_lip_beff)
        bot_K = Sec553.springStiffnessK(self.E, self.tcore, self.v, bot_b1, self.aa,
                                        bot_b1, True)
        bot_Is = Sec553.Is(bot_flg_be2, self.tcore, bot_lip_beff)
        bot_scrs = Sec553.calc_scrs(bot_K, bot_Is, self.E, bot_As)
        # Thickness reduction factor
        bot_xd = 1.0  # Tension on bottom flange
        bot_t_red = bot_xd * self.tcore
        # ==============================================================================================================
        # Effective width of the web
        # ==============================================================================================================
        # Stress ratio / Uniform compression on web
        ff = 1.0
        web_ksigma = Sec4.Table4_1_ksigma(ff)
        web_lamp = Sec4.lamp(self.aa, self.tcore, web_ksigma, scomed, True)
        web_rho = Sec4.internal_element(web_lamp, ff)
        web_beff = Sec4.Table4_1_beff(ff, self.aa, web_rho)[0]
        web_be1 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[1]
        web_be2 = Sec4.Table4_1_beff(ff, self.aa, web_rho)[2]
        h1 = web_be1
        h2 = web_be1
        # ==============================================================================================================
        # Effective section properties
        # ==============================================================================================================
        # Create a matrix contains the element data from bottom lip to top lip
        # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
        self.BendWeakWeb_elementData2 = np.array(
            [[1, self.bb, bot_lip_beff, self.bb, 0.0, bot_t_red],
             [2, self.bb, 0.0, self.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, self.tcore],
             [4, 0.0, 0.0, 0.0, h2, self.tcore],
             [5, 0.0, self.aa - h1, 0.0, self.aa, self.tcore],
             [6, 0.0, self.aa, top_flg_be1, self.aa, self.tcore],
             [7, self.bb - top_flg_be2, self.aa, self.bb, self.aa, top_t_red],
             [8, self.bb, self.aa, self.bb, self.aa - top_lip_beff, top_t_red]])

        # Results
        self.BendWeakWeb_Iyeff = intprop.calcProps(self.BendWeakWeb_elementData2)[4]
        self.BendWeakWeb_xgc = intprop.calcProps(self.BendWeakWeb_elementData2)[3]
        self.BendWeakWeb_xgct = self.bb - intprop.calcProps(self.BendWeakWeb_elementData2)[3]
        self.BendWeakWeb_Wyeff = intprop.calcProps(self.BendWeakWeb_elementData2)[4] / (
            max(self.BendWeakWeb_xgc, self.BendWeakWeb_xgct))

    def plot_effC_section(self):
        figure, axis = plt.subplots(2, 2)
        figure.suptitle("Effective Section")
        for i in self.Axial_elementData2:
            x = [i[1], i[3]]
            y = [i[2], i[4]]
            t = i[5]
            axis[0, 0].plot(x, y, linewidth=t * 2)
        axis[0, 0].set_title('Axial Compression')
        axis[0, 0].grid(color='green', linestyle='--', linewidth=0.5)
        axis[0, 0].plot(self.Axial_xgc, self.A / 2, color='gray', linestyle='solid', marker=6, markerfacecolor='blue',
                        markersize=9)
        axis[0, 0].text(self.Axial_xgc, self.A / 2 * 0.9, 'Effective Gravity Centre', style='italic',
                        bbox={'facecolor': 'blue', 'alpha': 0.5, 'pad': 1}, verticalalignment='top')
        axis[0, 0].plot(self.zgx, self.A / 2, color='gray', linestyle='solid', marker=7, markerfacecolor='red',
                        markersize=9)
        axis[0, 0].text(self.zgx, self.A / 2 * 1.2, 'Gross Gravity Centre', style='italic',
                        bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 1}, verticalalignment='top')
        axis[0, 0].text(self.zgx * 1.3, self.A / 2, f'eny : {self.Axial_dxgc:.2f} mm', style='italic',
                        bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 1}, verticalalignment='center')
        axis[0, 0].axis('equal')

        for i in self.BendStrong_elementData2:
            x = [i[1], i[3]]
            y = [i[2], i[4]]
            t = i[5]
            axis[0, 1].plot(x, y, linewidth=t * 2)
        axis[0, 1].set_title('Bending About Strong Axis')
        axis[0, 1].grid(color='green', linestyle='--', linewidth=0.5)
        axis[0, 1].axis('equal')

        for i in self.BendWeakLip_elementData2:
            x = [i[1], i[3]]
            y = [i[2], i[4]]
            t = i[5]
            axis[1, 0].plot(x, y, linewidth=t * 2)
        axis[1, 0].set_title(f'Bending About Weak Axis\nLips are Under Compression')
        axis[1, 0].grid(color='green', linestyle='--', linewidth=0.5)
        axis[1, 0].axis('equal')

        for i in self.BendWeakWeb_elementData2:
            x = [i[1], i[3]]
            y = [i[2], i[4]]
            t = i[5]
            axis[1, 1].plot(x, y, linewidth=t * 2)
        axis[1, 1].set_title(f'Bending About Weak Axis\nWeb is Under Compression')
        axis[1, 1].grid(color='green', linestyle='--', linewidth=0.5)
        axis[1, 1].axis('equal')
        plt.show()

    def plot_C_section(self):
        t = self.x
        s = self.y
        fig, ax = plt.subplots()
        ax.plot(t, s)
        ax.set(xlabel='mm', ylabel='mm', title=self.descp)
        ax.grid()
        plt.axis('equal')
        plt.show()





# Calculating the gross section properties
section = GrossProp(90.0, 45.0, 10.0, 1.2, 1.6, 350.0)
