import EffectiveSection.EN1993_1_5.Sec4 as Sec4
import EffectiveSection.EN1993_1_3.Sec5_5_3 as Sec553
import EffectiveSection.Modes.IntCalcProp as intprop
import Definitions.Definitions as defin
import numpy as np
import matplotlib.pyplot as plt


class bendWeakWeb:
    def __init__(self):
        # Variables for bending about strong axis
        self.BendWeakWeb_xgct = None
        self.BendWeakWeb_xgc = None
        self.BendWeakWeb_elementData2 = None
        self.BendWeakWeb_ygct = None  # Top flange extreme fiber to neutral axis
        self.BendWeakWeb_ygc = None  # Bottom flange extreme fiber to neutral axis
        self.BendWeakWeb_Iyeff = None  # Moment of inertia
        self.BendWeakWeb_Wyeff = None  # Section modulus
        self.calcs_BendingWeakWeb()

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

    def calcs_BendingWeakWeb(self):
        # Design Stress
        scomed = defin.steel.fy
        # ==============================================================================================================
        # Effective width of the top flange
        # ==============================================================================================================
        # Only tension part is taken into account for the flange (bb-zgx)
        top_flg_be1 = 0.0
        top_flg_be2 = (defin.section.bb - defin.gross.zgx)
        # ==============================================================================================================
        # Effective width of the bot flange
        # ==============================================================================================================
        # Only tension part is taken into account for the flange (bb-zgx)
        bot_flg_be1 = 0.0
        bot_flg_be2 = (defin.section.bb - defin.gross.zgx)
        # ==============================================================================================================
        # Effective width of the top edge fold
        # ==============================================================================================================
        # Lip is fully effective under tension.
        top_lip_stres_ratio = Sec4.stres_ratio(defin.section.cc, 0.0)  # bwhole is 0 for uniform stress.
        top_lip_ksigma = Sec553.ksig(defin.section.cc, defin.section.bb)
        top_lip_lamp = Sec4.lamp(defin.section.cc, defin.section.tcore, top_lip_ksigma, scomed, False)
        top_lip_rho = Sec4.outstand_element(top_lip_lamp)
        top_lip_beff = Sec4.Table4_2_beff(defin.section.cc, top_lip_rho, top_lip_stres_ratio)[0]
        top_lip_bc = Sec4.Table4_2_beff(defin.section.cc, top_lip_rho, top_lip_stres_ratio)[1]
        top_lip_bt = Sec4.Table4_2_beff(defin.section.cc, top_lip_rho, top_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        top_As = defin.section.tcore * (top_flg_be2 + top_lip_beff)
        # Spring stiffness of the edge stiffener
        top_b1 = Sec553.calc_b1(defin.section.bb, top_flg_be2, defin.section.tcore, top_lip_beff)
        top_K = Sec553.springStiffnessK(defin.steel.E, defin.section.tcore, defin.steel.v, top_b1, defin.section.aa,
                                        top_b1, True)
        top_Is = Sec553.Is(top_flg_be2, defin.section.tcore, top_lip_beff)
        top_scrs = Sec553.calc_scrs(top_K, top_Is, defin.steel.E, top_As)
        # Thickness reduction factor
        top_xd = Sec553.thk_reduction(scomed, top_scrs)
        top_t_red = top_xd * defin.section.tcore
        # ==============================================================================================================
        # Effective width of the bottom edge fold
        # ==============================================================================================================
        # Lip is fully effective under tension.
        bot_lip_stres_ratio = Sec4.stres_ratio(defin.section.cc, 0.0)  # bwhole is 0 for uniform stress.
        bot_lip_ksigma = Sec553.ksig(defin.section.cc, defin.section.bb)
        bot_lip_lamp = Sec4.lamp(defin.section.cc, defin.section.tcore, bot_lip_ksigma, scomed, False)
        bot_lip_rho = Sec4.outstand_element(bot_lip_lamp)
        bot_lip_beff = Sec4.Table4_2_beff(defin.section.cc, bot_lip_rho, bot_lip_stres_ratio)[0]
        bot_lip_bc = Sec4.Table4_2_beff(defin.section.cc, bot_lip_rho, bot_lip_stres_ratio)[1]
        bot_lip_bt = Sec4.Table4_2_beff(defin.section.cc, bot_lip_rho, bot_lip_stres_ratio)[2]
        # Effective area of the edge stiffener
        bot_As = defin.section.tcore * (bot_flg_be2 + bot_lip_beff)
        # Spring stiffness of the edge stiffener
        bot_b1 = Sec553.calc_b1(defin.section.bb, bot_flg_be2, defin.section.tcore, bot_lip_beff)
        bot_K = Sec553.springStiffnessK(defin.steel.E, defin.section.tcore, defin.steel.v, bot_b1, defin.section.aa,
                                        bot_b1, True)
        bot_Is = Sec553.Is(bot_flg_be2, defin.section.tcore, bot_lip_beff)
        bot_scrs = Sec553.calc_scrs(bot_K, bot_Is, defin.steel.E, bot_As)
        # Thickness reduction factor
        bot_xd = 1.0  # Tension on bottom flange
        bot_t_red = bot_xd * defin.section.tcore
        # ==============================================================================================================
        # Effective width of the web
        # ==============================================================================================================
        # Stress ratio / Uniform compression on web
        ff = 1.0
        web_ksigma = Sec4.Table4_1_ksigma(ff)
        web_lamp = Sec4.lamp(defin.section.aa, defin.section.tcore, web_ksigma, scomed, True)
        web_rho = Sec4.internal_element(web_lamp, ff)
        web_beff = Sec4.Table4_1_beff(ff, defin.section.aa, web_rho)[0]
        web_be1 = Sec4.Table4_1_beff(ff, defin.section.aa, web_rho)[1]
        web_be2 = Sec4.Table4_1_beff(ff, defin.section.aa, web_rho)[2]
        h1 = web_be1
        h2 = web_be1
        # ==============================================================================================================
        # Effective section properties
        # ==============================================================================================================
        # Create a matrix contains the element data from bottom lip to top lip
        # 0 id , 1 inodeX, 2 inodeY, 3 jnodeX, 4 JnodeY, 5 thickness
        self.BendWeakWeb_elementData2 = np.array(
            [[1, defin.section.bb, bot_lip_beff, defin.section.bb, 0.0, bot_t_red],
             [2, defin.section.bb, 0.0, defin.section.bb - bot_flg_be2, 0.0, bot_t_red],
             [3, bot_flg_be1, 0.0, 0.0, 0.0, defin.section.tcore],
             [4, 0.0, 0.0, 0.0, h2, defin.section.tcore],
             [5, 0.0, defin.section.aa - h1, 0.0, defin.section.aa, defin.section.tcore],
             [6, 0.0, defin.section.aa, top_flg_be1, defin.section.aa, defin.section.tcore],
             [7, defin.section.bb - top_flg_be2, defin.section.aa, defin.section.bb, defin.section.aa, top_t_red],
             [8, defin.section.bb, defin.section.aa, defin.section.bb, defin.section.aa - top_lip_beff, top_t_red]])

        for i in self.BendWeakWeb_elementData2:
            x = [i[1], i[3]]
            y = [i[2], i[4]]
            t = i[5]
            plt.plot(x, y, linewidth=t*2)
        plt.title('Bending About Weak Axis Web is Under Compression | Effective Section')
        plt.axis('equal')
        plt.show()

        # Results
        self.BendWeakWeb_Iyeff = intprop.calcProps(self.BendWeakWeb_elementData2)[4]
        self.BendWeakWeb_xgc = intprop.calcProps(self.BendWeakWeb_elementData2)[3]
        self.BendWeakWeb_xgct = defin.section.bb - intprop.calcProps(self.BendWeakWeb_elementData2)[3]
        self.BendWeakWeb_Wyeff = intprop.calcProps(self.BendWeakWeb_elementData2)[4] / (
            max(self.BendWeakWeb_xgc, self.BendWeakWeb_xgct))


Bending_weak_web = bendWeakWeb()
