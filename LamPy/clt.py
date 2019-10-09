import numpy as np
import math as mt
import matplotlib.pyplot as plt
from tikzplotlib import save as tikz_save
from matplotlib.patches import Rectangle, Polygon, Patch
import pandas as pd
from matplotlib import cm

np.seterr(divide='ignore', invalid='ignore')


class laminate(object):
    def __init__(self, Interactive=False, PrintAll=True, Plot=True,
                 nPly=3, thkPly=[1, 1, 1], theta=np.array([-45, -30, 0]),
                 E1=141000, E2=8500, G12=3750, nu12=0.25, phi=0.6,
                 EFiber=230000, EMatrix=3500, nu12Fiber=0.2, nu12Matrix=0.33,
                 symmetric=True, RhoFiber=1.8, RhoMatrix=1.2,  Nx=100, Ny=0,
                 Nxy=0, Mx=0, My=0, Mxy=0, T1=1950, T2=48, C1=1480, C2=200,
                 S12=55, S13=55, S23=55, LaminateName=""):
        self.Interactive = Interactive
        self.PrintAll = PrintAll
        self.Plot = Plot
        self.symmetric = symmetric
        self.nPly = nPly
        self.theta = theta
        self.EFiber = EFiber
        self.EMatrix = EMatrix
        self.nu12Fiber = nu12Fiber
        self.nu12Matrix = nu12Matrix
        self.E1 = E1
        self.E2 = E2
        self.G12 = G12
        self.nu12 = nu12
        self.phi = phi         # fiber volume percentage
        self.thkPly = thkPly
        self.Nx = Nx
        self.Ny = Ny
        self.Nxy = Nxy
        self.Mx = Mx
        self.My = My
        self.Mxy = Mxy
        self.T1 = T1
        self.T2 = T2
        self.C1 = C1
        self.C2 = C2
        self.S12 = S12
        self.S13 = S12
        self.S23 = S23
        self.RhoFiber = RhoFiber
        self.RhoMatrix = RhoMatrix
        self.LaminateName = LaminateName
        self.unitAngles = "deg"

    def calcRuleOfMixtures(self):
        self.E1 = self.EFiber*self.phi+self.EMatrix*(1-self.phi)
        self.E2 = (self.phi/self.EFiber + (1-self.phi)/self.EMatrix)**-1
        self.G12Matrix = self.EMatrix/(2*(1+self.nu12Matrix))
        self.G12Fiber = self.EFiber/(2*(1+self.nu12Fiber))
        self.G12 = (self.phi/self.G12Fiber +
                    (1-self.phi)/self.G12Matrix)**-1
#        self.G12 = self.G12Matrix*self.G12Fiber/((self.G12Matrix*self.phi)+
#                                                 (1-self.phi)*self.G12Fiber)
        self.nu12 = self.nu12Fiber*self.phi+(1-self.phi)*self.nu12Matrix

    def SymmetricLaminate(self):
        if self.symmetric:
            ThkSym = np.concatenate((self.thkPly, np.flip(self.thkPly)),
                                    axis=None)
            ThetaSym = np.concatenate((self.theta, np.flip(-self.theta)),
                                      axis=None)
            PlySym = 2*self.nPly
            self.thkPly = ThkSym
            self.theta = ThetaSym
            self.nPly = PlySym

    def calcPlyPositions(self):
        self.thkLam = sum(self.thkPly)
        self.z = [[]]*self.nPly
        self.h = [[]]*self.nPly
        for i in range(self.nPly):
            self.z[i] = np.array([self.thkLam/2-sum(self.thkPly[0:i]),
                                  -self.thkLam/2+sum(self.thkPly[i+1:])])
            self.h[i] = self.z[i][np.argmax((np.abs([self.z[i][0],
                                                     self.z[i][1]])))]

    def calcPoisson21(self):
        self.nu21 = (self.nu12*self.E2)/self.E1

    def deg2rad(self):
        if self.unitAngles == "deg":
            self.theta = np.deg2rad(self.theta)
#            self.unitAngles = "rad"

    def calcABD(self):
        Q11 = self.E1/(1-(self.E2*self.nu12**2)/self.E1)
        Q22 = self.E2/self.E1*Q11
        Q12 = self.nu12*Q22
        Q66 = self.G12
        Q = np.array([[Q11, Q12,   0],
                      [Q12, Q22,   0],
                      [  0,   0, Q66]])
        self.A = np.zeros((3, 3))
        self.B = np.zeros((3, 3))
        self.D = np.zeros((3, 3))
        self.Qbar = [[]]*self.nPly
        self.T = [[]]*self.nPly
        self.Tinv = [[]]*self.nPly
        for i in range(self.nPly):
            c = np.cos((self.theta[i]))
            s = np.sin((self.theta[i]))
            cs = c*s
            cc = c**2
            ss = s**2
            self.T[i] = np.array([[ cc, ss,  2*cs],
                                  [ ss, cc, -2*cs],
                                  [-cs, cs,  cc-ss]])
            self.Tinv[i] = np.linalg.inv(self.T[i])
            self.Qbar[i] = self.Tinv[i] @ Q @ self.Tinv[i].T
            self.A += self.Qbar[i]*(self.z[i][0] - self.z[i][1])
            self.B += 1/2*self.Qbar[i]*(self.z[i][0]**2 - self.z[i][1]**2)
            self.D += 1/3*self.Qbar[i]*(self.z[i][0]**3 - self.z[i][1]**3)
        self.ABD = np.block([[self.A, self.B],
                             [self.B, self.D]])

    def cleanABD(self):
        self.ABD[abs(self.ABD) < 1e-3] = 0.0

    def calcInverseABD(self):
        self.ABDInv = np.linalg.inv(self.ABD)

    def calcLaminateParameters(self):
        self.ELx = 1/(self.ABDInv[0, 0]*self.thkLam)
        self.ELy = 1/(self.ABDInv[1, 1]*self.thkLam)
        self.GLxy = 1/(self.ABDInv[2, 2]*self.thkLam)
        self.nuLxy = self.ABDInv[0, 1]/self.ABDInv[0, 0]
        self.nuLyx = self.ABDInv[0, 1]/self.ABDInv[1, 1]

    def calcStrain(self):
        loadVec = np.array([[self.Nx, self.Ny, self.Nxy, self.Mx, self.My,
                             self.Mxy]]).T
        self.strainVec = self.ABDInv@loadVec
        #self.strainLVec = [[]]*self.nPly
        self.strainLTopVec = [[]]*self.nPly
        self.strainLBotVec = [[]]*self.nPly
        self.strainPlyTopVec = [[]]*self.nPly
        self.strainPlyBotVec = [[]]*self.nPly
        for i in range(self.nPly):
            if self.z[i][0] > 0:
                self.strainLTopVec[i] = (self.strainVec[:3] +
                                         self.strainVec[3:]*(self.z[i][0]))
                self.strainLBotVec[i] = (self.strainVec[:3] +
                                         self.strainVec[3:]*(self.z[i][1]))
            else:
                self.strainLTopVec[i] = (self.strainVec[:3] +
                                         self.strainVec[3:]*(self.z[i][1]))
                self.strainLBotVec[i] = (self.strainVec[:3] +
                                         self.strainVec[3:]*(self.z[i][0]))
            self.strainPlyTopVec[i] = (np.linalg.inv(self.T[i].T) @
                                       self.strainLTopVec[i])
            self.strainPlyBotVec[i] = (np.linalg.inv(self.T[i].T) @
                                       self.strainLBotVec[i])

    def calcStress(self):
        self.stressLTopVec = [[]]*self.nPly
        self.stressLBotVec = [[]]*self.nPly
        self.stressPlyTopVec = [[]]*self.nPly
        self.stressPlyBotVec = [[]]*self.nPly
        for i in range(self.nPly):
            self.stressLTopVec[i] = self.Qbar[i] @ self.strainLTopVec[i]
            self.stressLBotVec[i] = self.Qbar[i] @ self.strainLBotVec[i]
            self.stressPlyTopVec[i] = self.T[i] @ self.stressLTopVec[i]
            self.stressPlyBotVec[i] = self.T[i] @ self.stressLTopVec[i]

    def calcMisesStress(self):
        A = np.array([[ 1.0, -0.5, 0.0],
                      [-0.5,  1.0, 0.0],
                      [ 0.0,  0.0, 3.0]])
        self.stressMisesTop = np.zeros((self.nPly,))
        self.stressMisesBot = np.zeros((self.nPly,))
        for i in range(self.nPly):
            self.stressMisesBot[i] = np.sqrt(self.stressLBotVec[i].T @ A @
                                             self.stressLBotVec[i])
            self.stressMisesTop[i] = np.sqrt(self.stressLTopVec[i].T @ A @
                                             self.stressLTopVec[i])

    def calcFailureTsaiWu(self):
        self.FailureTsaiWu = []
        self.ReserveTsaiWu = []
        F11 = 1/(self.T1*self.C1)
        F22 = 1/(self.T2*self.C2)
        F1 = (1/self.T1) - (1/self.C1)
        F2 = (1/self.T2) - (1/self.C2)
        F12 = -0.5/np.sqrt(self.T1*self.T2*self.C1*self.C2)
        F33 = 1/(self.S12**2)
        self.FailureModeTsaiWu = [[]]*self.nPly
        self.ReserveTsaiWu = [[]]*self.nPly
        for i in range(self.nPly):
            Mode = [[]]*3
            aTop = (F11*self.stressPlyTopVec[i][0]**2 +
                    F22*self.stressPlyTopVec[i][1]**2 +
                    F33*self.stressPlyTopVec[i][2]**2 +
                    2*F12*self.stressPlyTopVec[i][0] *
                    self.stressPlyTopVec[i][1])[0]
            aBot = (F11*self.stressPlyBotVec[i][0]**2 +
                    F22*self.stressPlyBotVec[i][1]**2 +
                    F33*self.stressPlyBotVec[i][2]**2 +
                    2*F12*self.stressPlyBotVec[i][0] *
                    self.stressPlyBotVec[i][1])[0]
            bTop = (F1*self.stressPlyTopVec[i][0] +
                    F2*self.stressPlyTopVec[i][1])[0]
            bBot = (F1*self.stressPlyBotVec[i][0] +
                    F2*self.stressPlyBotVec[i][1])[0]
            a = max(aTop, aBot)
            b = max(bTop, bBot)
            self.ReserveTsaiWu[i] = (-b + (b**2 + 4*a)**(1/2))/(2*a)
            if self.ReserveTsaiWu[i] < 1.0:
                Modes = ["Fiber", "Matrix", "Shear"]
                Mode[0] = max(abs(F1*self.stressPlyTopVec[i][0] +
                                  F11*self.stressPlyTopVec[i][0]**2),
                              abs(F1*self.stressPlyBotVec[i][0] +
                                  F11*self.stressPlyBotVec[i][0]**2))
                Mode[1] = max(abs(F2*self.stressPlyTopVec[i][1] +
                                  F22*self.stressPlyTopVec[i][1]**2),
                              abs(F2*self.stressPlyBotVec[i][1] +
                                  F22*self.stressPlyBotVec[i][1]**2))
                Mode[2] = max(abs(F33*self.stressPlyTopVec[i][2]**2),
                              abs(F33*self.stressPlyBotVec[i][2]**2))
                self.FailureModeTsaiWu[i] = Modes[np.argmax(Mode)]
            else:
                self.FailureModeTsaiWu[i] = "None"

    def calcMass(self):
        self.LaminateVolume = (np.sum(self.thkPly)/1000)  # fibers volume in m^3 for one m^2 of laminate
        self.LaminateMass = self.LaminateVolume*(self.RhoFiber*self.phi +
                                                 self.RhoMatrix*(1-self.phi))
        #Laminate mass in tons, referred to one m^2

    def plotStackLayup(self, savePng=True, saveSvg=True, saveTex=True,
                       plotShow=True, fontSize=12,
                       fontName="Tex Gyre Pagella", ColorMap=0):
        #if fontName == "Palatino":
        #    fontName = "Tex Gyre Pagella"
        colorList = defineColorMap(ColorMap)
        fig = plt.figure(figsize=(2+int(self.nPly-1)/2, 2+int(self.nPly-1)/2))
        ax = fig.add_subplot()

        def plot_hatches(ax, angle, ii, delta, offset=0.0):
            nPoints = 501
            xL = -2
            xU = 2
            x = np.linspace(xL, xU, nPoints)
            for c in np.arange(xL, xU, offset):
                yprime = mt.cos(angle) * c + mt.sin(angle) * x + delta
                xprime = -mt.sin(angle) * c + mt.cos(angle) * x + delta
                yprime[yprime > (0.999+delta)] = np.nan
                xprime[xprime > (0.999+delta)] = np.nan
                yprime[yprime < (0.001+delta)] = np.nan
                yprime[yprime < (0.001+delta)] = np.nan
                xprime[xprime < (0.001+delta)] = np.nan
                xprimeTemp = xprime[np.logical_not(np.isnan(xprime))]
                yprimeTemp = yprime[np.logical_not(np.isnan(xprime))]
                xprime = xprimeTemp[np.logical_not(np.isnan(yprimeTemp))]
                yprime = yprimeTemp[np.logical_not(np.isnan(yprimeTemp))]
                ax.plot(xprime, yprime, color="k", zorder=ii)
            return ax
        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly-1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = i
            index[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = -1
        for ii in range(self.nPly):
            el = ii*0.25
            plot_hatches(ax, angle[ii], ii, el, offset=0.05)
            ax.add_patch(Rectangle((0+el, 0+el), 1, 1, fill="r", alpha=0.85,
                                   color=colorList[int(colorIndex[ii])],
                                   zorder=ii))
            plt.text(1.05+el, 0.05+el,
                     str(round(formatNumber(round(angle[ii], 4)), 2))+" deg",
                     fontsize=fontSize)
        ax.set_ylim(0, 1+el)
        ax.set_xlim(0, 1+el)
        plt.axis('off')
        #plt.rcParams["font.family"] = fontName
        if self.LaminateName == "":
            Name = ""
        else:
            Name = self.LaminateName+"_"
        if savePng:
            plt.savefig(Name + "StackLayup.png", format="png",
                        bbox_inches="tight")
        if saveSvg:
            plt.savefig(Name + "StackLayup.svg", format="svg",
                        bbox_inches="tight")
        if saveTex:  # Currently not working correctly, patches not shown.
            tikz_save(Name + "StackLayup.tex", show_info=False, strict=False,
                      extra_axis_parameters={"axis lines*=left"})
        if plotShow:
            plt.show()

    def plotPlyStrainStress(self, Show=True, SaveTex=True, SavePng=True,
                            SaveSvg=True, ColorMap=0, Legend=True):

        colorList = defineColorMap(ColorMap)
        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly-1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = i
            index[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = -1
        legendElements = []
        angleList = []
        angleOrder = np.argsort(angle)
        for i in range(len(angle)):
            j = angleOrder[i]
            val = abs(self.theta[j])
            if val not in angleList:
                angleList.append(val)
                lab = (str(round(formatNumber(round(abs(angle[j]), 4)), 2)) +
                       " [deg]")
                legendElements.append(Patch(facecolor=colorList[int(colorIndex[j])],
                                            edgecolor=colorList[int(colorIndex[j])],
                                            label=lab))
        fig = plt.figure(figsize=(7, 3))
        ax = fig.add_subplot()
        width = 1.5
        xMin = -width/2
        Shift = 2
        xShift = []
        x = [[]]*self.nPly
        for i in range(self.nPly):
            ax.add_patch(Rectangle((xMin, self.z[i][0]), width,
                                   -self.thkPly[i], fill="r", alpha=0.85,
                                   color=colorList[int(colorIndex[i])],
                                   zorder=i))
        for ii in range(3):
            xShift.append(1*Shift+2*Shift*ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (np.array([min(self.strainLTopVec[i][ii][0], 0),
                                      max(self.strainLTopVec[i][ii][0], 0),
                                      max(self.strainLBotVec[i][ii][0], 0),
                                      min(self.strainLBotVec[i][ii][0], 0)]) /
                            np.max(np.max(self.strainLTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
                else:
                    x[i] = (np.array([max(self.strainLBotVec[i][ii][0], 0),
                                      min(self.strainLBotVec[i][ii][0], 0),
                                      min(self.strainLTopVec[i][ii][0], 0),
                                      max(self.strainLTopVec[i][ii][0], 0)]) /
                            np.max(np.max(self.strainLTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
            xShift.append(2*Shift+2*Shift*ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (np.array([min(self.stressLTopVec[i][ii][0], 0),
                                      max(self.stressLTopVec[i][ii][0], 0),
                                      max(self.stressLBotVec[i][ii][0], 0),
                                      min(self.stressLBotVec[i][ii][0], 0)]) /
                            np.max(np.max(self.stressLTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
                else:
                    x[i] = (np.array([min(self.stressLBotVec[i][ii][0], 0),
                                      max(self.stressLBotVec[i][ii][0], 0),
                                      max(self.stressLTopVec[i][ii][0], 0),
                                      min(self.stressLTopVec[i][ii][0], 0)]) /
                            np.max(np.max(self.stressLTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
        plt.ylabel("laminate height [mm]")
        ax.set_ylim([np.min(self.z), np.max(self.z)])
        #plt.xticks(xShift, ('0\n strain x', '0\n stress x', '0\n strain y', '0\n stress y',
        #                    '0\nstrain xy', '0\n stress xy'))
        #plt.xticks(xShift, ['0', '0', '0','0', '0', '0'])
        # HERE AN ATTEMPT TO MAKE LABELS WORK IN TIKZ...
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")
        ax.set_xticks(xShift)
        ax.set_xticklabels(('0', '0', '0', '0', '0', '0'))
        plt.xlabel("values normalized based on maximum")
        ax2 = ax.twiny()
        ax2.xaxis.set_tick_params(labeltop='on')
        #ax2.xaxis.set_ticks_position("top")
        #ax2.xaxis.set_label_position("top")
        #ax2.spines["top"].set_position(("axes", -0.25))
        #ax2.set_frame_on(False)
        #ax2.patch.set_visible(False)
        #ax2.spines["top"].set_visible(False)
        ax2.set_xticks(np.insert(xShift, 0, 0))
        ax2.set_xticklabels(('stackup', 'strain x', 'stress x', 'strain y',
                             'stress y', 'strain xy', 'stress xy'))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.tick_params(axis=u'both', which=u'both', length=0)
        ax2.set_xlim([xMin, 13])
        ax.set_xlim([xMin, 13])
        if Legend:
            plt.legend(handles=legendElements[::-1], loc='center left',
                       bbox_to_anchor=(1.05, 0.5), frameon=False)
        if SaveTex:
            tikz_save('LaminateStrainStress.tex', show_info=False,
                      strict=False)  # ,figureheight='\\figureheight', figurewidth='\\figurewidth', extra_axis_parameters={"axis lines*=left"})
        if SavePng:
            plt.savefig("LaminateStrainStress.png", format="png")
        if SaveSvg:
            plt.savefig("LaminateStrainStress.svg", format="svg")
        if Show:
            plt.show()

    def plotPlyStrainStressPly(self, Show=True, SaveTex=True, SavePng=True,
                               SaveSvg=True, ColorMap=0, Legend=True):
        colorList = defineColorMap(ColorMap)
        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly-1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = i
            index[(np.abs(np.abs(angle)-np.abs(angle[i])) <= 0.001).tolist()] = -1
        legendElements = []
        angleList = []
        angleOrder = np.argsort(angle)
        for i in range(len(angle)):
            j = angleOrder[i]
            val = abs(self.theta[j])
            if val not in angleList:
                angleList.append(val)
                lab = (str(round(formatNumber(round(abs(angle[j]), 4)), 2)) +
                       " [deg]")
                legendElements.append(Patch(facecolor=colorList[int(colorIndex[j])],
                                            edgecolor=colorList[int(colorIndex[j])],
                                            label=lab))
        fig = plt.figure(figsize=(7, 3))
        ax = fig.add_subplot()
        width = 1.5
        xMin = -width/2
        Shift = 2
        xShift = []
        x = [[]]*self.nPly
        for i in range(self.nPly):
            ax.add_patch(Rectangle((xMin, self.z[i][0]), width,
                                   -self.thkPly[i], fill="r", alpha=0.85,
                                   color=colorList[int(colorIndex[i])],
                                   zorder=i))
        for ii in range(3):
            xShift.append(1*Shift+2*Shift*ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (np.array([min(self.strainPlyTopVec[i][ii][0], 0),
                                      max(self.strainPlyTopVec[i][ii][0], 0),
                                      max(self.strainPlyBotVec[i][ii][0], 0),
                                      min(self.strainPlyBotVec[i][ii][0], 0)])/
                            np.max(np.max(self.strainPlyTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
                else:
                    x[i] = (np.array([max(self.strainPlyBotVec[i][ii][0], 0),
                                      min(self.strainPlyBotVec[i][ii][0], 0),
                                      min(self.strainPlyTopVec[i][ii][0], 0),
                                      max(self.strainPlyTopVec[i][ii][0], 0)])/
                            np.max(np.max(self.strainPlyTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
            xShift.append(2*Shift+2*Shift*ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (np.array([min(self.stressPlyTopVec[i][ii][0], 0),
                                      max(self.stressPlyTopVec[i][ii][0], 0),
                                      max(self.stressPlyBotVec[i][ii][0], 0),
                                      min(self.stressPlyBotVec[i][ii][0], 0)])/
                            np.max(np.max(self.stressPlyTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
                else:
                    x[i] = (np.array([min(self.stressPlyBotVec[i][ii][0], 0),
                                      max(self.stressPlyBotVec[i][ii][0], 0),
                                      max(self.stressPlyTopVec[i][ii][0], 0),
                                      min(self.stressPlyTopVec[i][ii][0], 0)])/
                            np.max(np.max(self.stressPlyTopVec))+xShift[-1])
                    y = np.array([self.z[i][0], self.z[i][0], self.z[i][1],
                                  self.z[i][1]])
                    ax.add_patch(Polygon(xy=list(zip(x[i], y)), fill="r",
                                         alpha=0.85,
                                         color=colorList[int(colorIndex[i])],
                                         zorder=i))
        plt.ylabel("laminate height [mm]")
        ax.set_ylim([np.min(self.z), np.max(self.z)])
        #plt.xticks(xShift, ('0\n strain x', '0\n stress x', '0\n strain y', '0\n stress y',
        #                    '0\nstrain xy', '0\n stress xy'))
        #plt.xticks(xShift, ['0', '0', '0','0', '0', '0'])
        # HERE AN ATTEMPT TO MAKE LABELS WORK IN TIKZ...
        ax.xaxis.set_ticks_position("bottom")
        ax.xaxis.set_label_position("bottom")
        ax.set_xticks(xShift)
        ax.set_xticklabels(('0', '0', '0', '0', '0', '0'))
        plt.xlabel("values normalized based on maximum")
        ax2 = ax.twiny()
        ax2.xaxis.set_tick_params(labeltop='on')
        #ax2.xaxis.set_ticks_position("top")
        #ax2.xaxis.set_label_position("top")
        #ax2.spines["top"].set_position(("axes", -0.25))
        #ax2.set_frame_on(False)
        #ax2.patch.set_visible(False)
        #ax2.spines["top"].set_visible(False)
        ax2.set_xticks(np.insert(xShift, 0, 0))
        ax2.set_xticklabels(('stackup', 'strain 1', 'stress 1', 'strain 2',
                             'stress 2', 'strain 12', 'stress 12'))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax2.tick_params(axis=u'both', which=u'both', length=0)
        ax2.set_xlim([xMin, 13])
        ax.set_xlim([xMin, 13])
        if Legend:
            plt.legend(handles=legendElements[::-1], loc='center left',
                       bbox_to_anchor=(1.05, 0.5), frameon=False)
        if SaveTex:
            tikz_save('PlyStrainStress.tex', show_info=False, strict=False)
                 # figureheight='\\figureheight', figurewidth='\\figurewidth',
                 # extra_axis_parameters={"axis lines*=left"})
        if SavePng:
            plt.savefig("PlyStrainStress.png", format="png")
        if SaveSvg:
            plt.savefig("PlyStrainStress.svg", format="svg")
        if Show:
            plt.show()

    def writeTablePly(self, writeTex=True, writeCsv=True, printData=True,
                      DecimalPlaces=4, DecimalSymbol=".",
                      LaminateCooSys=True, PlyCooSys=True,
                      PlyTop=True, PlyBottom=True):
        FileName = "PlyTable"
        TestDataTop = {"ply": np.arange(1, self.nPly+1),
                       "fiber orientation": np.rad2deg(self.theta),
                       "strain x": [i[0][0] for i in self.strainLTopVec],
                       "stress x": [i[0][0] for i in self.stressLTopVec],
                       "strain y": [i[1][0] for i in self.strainLTopVec],
                       "stress y": [i[1][0] for i in self.stressLTopVec],
                       "strain xy": [i[2][0] for i in self.strainLTopVec],
                       "stress xy": [i[2][0] for i in self.stressLTopVec],
                       "strain 1": [i[0][0] for i in self.strainPlyTopVec],
                       "stress 1": [i[0][0] for i in self.stressPlyTopVec],
                       "strain 2": [i[1][0] for i in self.strainPlyTopVec],
                       "stress 2": [i[1][0] for i in self.stressPlyTopVec],
                       "strain 12": [i[2][0] for i in self.strainPlyTopVec],
                       "stress 12": [i[2][0] for i in self.stressPlyTopVec],
                       "Tsai-Wu reserve": self.ReserveTsaiWu}
        TestDataBot = {"ply": np.arange(1, self.nPly+1),
                       "fiber orientation": np.rad2deg(self.theta),
                       "strain x": [i[0][0] for i in self.strainLBotVec],
                       "stress x": [i[0][0] for i in self.stressLBotVec],
                       "strain y": [i[1][0] for i in self.strainLBotVec],
                       "stress y": [i[1][0] for i in self.stressLBotVec],
                       "strain xy": [i[2][0] for i in self.strainLBotVec],
                       "stress xy": [i[2][0] for i in self.stressLBotVec],
                       "strain 1": [i[0][0] for i in self.strainPlyBotVec],
                       "stress 1": [i[0][0] for i in self.stressPlyBotVec],
                       "strain 2": [i[1][0] for i in self.strainPlyBotVec],
                       "stress 2": [i[1][0] for i in self.stressPlyBotVec],
                       "strain 12": [i[2][0] for i in self.strainPlyBotVec],
                       "stress 12": [i[2][0] for i in self.stressPlyBotVec],
                       "Tsai-Wu reserve": self.ReserveTsaiWu}
        ColumnLam = ["ply", "fiber orientation", "strain x", "stress x",
                     "strain y",  "stress y", "strain xy",  "stress xy",
                     "Tsai-Wu reserve"]
        ColumnPly = ["ply", "fiber orientation", "strain 1", "stress 1",
                     "strain 2", "stress 2", "strain 12", "stress 12",
                     "Tsai-Wu reserve"]
        LaminateDataFrameTop = pd.DataFrame(TestDataTop, columns=ColumnLam)
        LaminateDataFrameBot = pd.DataFrame(TestDataBot, columns=ColumnLam)
        PlyDataFrameTop = pd.DataFrame(TestDataTop, columns=ColumnPly)
        PlyDataFrameBot = pd.DataFrame(TestDataBot, columns=ColumnPly)
        LaminateDataFrameTop = LaminateDataFrameTop.round(DecimalPlaces)
        LaminateDataFrameBot = LaminateDataFrameBot.round(DecimalPlaces)
        PlyDataFrameTop = PlyDataFrameTop.round(DecimalPlaces)
        PlyDataFrameBot = PlyDataFrameBot.round(DecimalPlaces)
        if writeTex:
            if LaminateCooSys:
                headerTex = ["ply", "fiber orientation $\\theta$ [deg]",
                             "strain x $\\varepsilon_x$ [-]",
                             "stress x $\\sigma_x$ [MPa]",
                             "strain y $\\varepsilon_y$ [-]",
                             "stress y $\\sigma_y$ [MPa]",
                             "strain xy $\\varepsilon_{xy}$ [-]",
                             "stress xy $\\sigma_{xy}$ [MPa]",
                             "Tsai-Wu reserve $R_\\mathrm{TW}$ [-]"]
                if PlyTop:
                    LaminateDataFrameTop.to_latex(FileName+"LaminateTop.tex",
                                                  index=None,
                                                  decimal=DecimalSymbol,
                                                  header=headerTex,
                                                  escape=False)
                    if writeCsv:
                        LaminateDataFrameTop.to_excel(FileName +
                                                      "LaminateTop.xlsx",
                                                      index=None, header=True)
                    if printData:
                        print(self.LaminateName,
                              'Laminate top reference system characteristics')
                        print(LaminateDataFrameTop.to_string(index=False))

                if PlyBottom:
                    LaminateDataFrameBot.to_latex(FileName+"LaminateBot.tex",
                                                  index=None,
                                                  decimal=DecimalSymbol,
                                                  header=headerTex,
                                                  escape=False)
                    if writeCsv:
                        LaminateDataFrameBot.to_excel(FileName +
                                                      "LaminateBot.xlsx",
                                                      index=None, header=True)
                    if printData:
                        print(self.LaminateName,
                              'Laminate bottom reference system characteristics')
                        print(LaminateDataFrameBot.to_string(index=False))
            if PlyCooSys:
                headerTex = ["ply", "fiber orientation $\\theta$ [deg]",
                             "strain 1 $\\varepsilon_1$ [-]",
                             "stress 1 $\\sigma_1$ [MPa]",
                             "strain 2 $\\varepsilon_2$ [-]",
                             "stress 2 $\\sigma_2$ [MPa]",
                             "strain 12 $\\varepsilon_{12}$ [-]",
                             "stress 12 $\\sigma_{12}$ [MPa]",
                             "Tsai-Wu reserve $R_\\mathrm{TW}$ [-]"]
                if PlyTop:
                    PlyDataFrameTop.to_latex(FileName+"PlyTop.tex", index=None,
                                             decimal=DecimalSymbol,
                                             header=headerTex, escape=False)
                    if writeCsv:
                        PlyDataFrameTop.to_excel(FileName+"PlyTop.xlsx",
                                                 index=None, header=True)
                    if printData:
                        print(self.LaminateName,
                              'Ply top reference system characteristics')
                        print(PlyDataFrameTop.to_string(index=False))
                if PlyBottom:
                    PlyDataFrameBot.to_latex(FileName+"PlyBot.tex", index=None,
                                             decimal=DecimalSymbol,
                                             header=headerTex, escape=False)
                    if writeCsv:
                        PlyDataFrameBot.to_excel(FileName+"PlyBot.xlsx",
                                                 index=None, header=True)
                    if printData:
                        print(self.LaminateName,
                              'Ply bottom reference system characteristics')
                        print(PlyDataFrameBot.to_string(index=False))


def formatNumber(value):
    if value % 1 == 0:
        value = int(value)
    return(value)


def plotParameterStudyAngle(angles, TW_Res, VMS, Ep, strain, Show=True,
                            SaveTex=True, SavePng=True, SaveSvg=True, Name="",
                            buffer=0.1, Grid=False, plotSize=(7, 5),
                            ColorTW="r", SymbolTW=".", LineTW="-",
                            ColorVM="b", SymbolVM=".", LineVM="-",
                            ColorStrainFiber="g", SymbolStrainFiber=".",
                            LineStrainFiber="-", ColorStrainx="g",
                            SymbolStrainx=".", LineStrainx="-", xBuffer=0):
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    plt.plot(angles, TW_Res, ColorTW+SymbolTW+LineTW)
    plt.ylabel('Tsai-Wu reserve factor $R_{TW}$ [-]')
    plt.xlim([min(angles)-xBuffer, max(angles)+xBuffer])
    plt.ylim([min(TW_Res)-max(np.abs(TW_Res))*buffer,
              max(TW_Res)+max(np.abs(TW_Res))*buffer])
    yticks = np.arange(0, max(TW_Res)+10, step=10)
    yticks[0] = 1
    plt.yticks(yticks)
    xticks = np.arange(min(angles), max(angles)+1, step=15)
    plt.xticks(xticks)
    plt.xlabel("fiber orientation angle $\\theta$ [deg]")
    if SaveTex:
        tikz_save(Name+'_TsaiWuReserve.tex', show_info=False, strict=False,
                  figureheight='\\figureheight', figurewidth='\\figurewidth',
                  extra_axis_parameters={"axis lines*=left"})
    if SavePng:
        plt.savefig(Name+"_TsaiWuReserve.png", format="png")
    if SaveSvg:
        plt.savefig(Name+"_TsaiWuReserve.svg", format="svg")
    if Show:
        plt.show()
    fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')
    plt.plot(angles, VMS, ColorVM+SymbolVM+LineVM)
    plt.ylabel('equivalent stress after Mises $\\sigma_M$ [MPa] ')
    plt.xlim([min(angles)-xBuffer, max(angles)+xBuffer])
    plt.ylim([min(VMS)-max(np.abs(VMS))*buffer,
              max(VMS)+max(np.abs(VMS))*buffer])
    xticks = np.arange(min(angles), max(angles)+1, step=15)
    plt.xticks(xticks)
    plt.xlabel("fiber orientation angle $\\theta$ [deg]")
    plt.tight_layout()
    if SaveTex:
        tikz_save(Name+'_StressEqv.tex', show_info=False, strict=False,
                  figureheight='\\figureheight', figurewidth='\\figurewidth',
                  extra_axis_parameters={"axis lines*=left"})
    if SavePng:
        plt.savefig(Name+"_StressMises.png", format="png")
    if SaveSvg:
        plt.savefig(Name+"_StressMises.svg", format="svg")
    if Show:
        plt.show()
    fig, ax3 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.yaxis.set_ticks_position('left')
    ax3.xaxis.set_ticks_position('bottom')
    plt.plot(angles, Ep, ColorStrainx+SymbolStrainx+LineStrainx)
    plt.ylabel('strain in $x$ [-]')
    plt.xlim([min(angles)-xBuffer, max(angles)+xBuffer])
    plt.ylim([min(Ep)-max(np.abs(Ep))*buffer,
              max(Ep)+max(np.abs(Ep))*buffer])
    xticks = np.arange(min(angles), max(angles)+1, step=15)
    plt.xticks(xticks)
    plt.xlabel("fiber orientation angle $\\theta$ [deg]")
    plt.tight_layout()
    if SaveTex:
        tikz_save(Name+'_Strain_x.tex', show_info=False, strict=False,
                  figureheight='\\figureheight', figurewidth='\\figurewidth',
                  extra_axis_parameters={"axis lines*=left"})
    if SavePng:
        plt.savefig(Name+"_Strain_x.png", format="png")
    if SaveSvg:
        plt.savefig(Name+"_Strain_x.svg", format="svg")
    if Show:
        plt.show()
    fig, ax4 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.yaxis.set_ticks_position('left')
    ax4.xaxis.set_ticks_position('bottom')
    plt.plot(angles, strain,
             ColorStrainFiber+SymbolStrainFiber+LineStrainFiber)
    plt.ylabel('strain in $\\parallel$ [-]')
    plt.xlim([min(angles)-xBuffer, max(angles)+xBuffer])
    plt.ylim([min(strain)-max(np.abs(strain))*buffer,
              max(strain)+max(np.abs(strain))*buffer])
    xticks = np.arange(min(angles), max(angles)+1, step=15)
    plt.xticks(xticks)
    plt.xlabel("fiber orientation angle $\\theta$ [deg]")
    plt.tight_layout()
    if SaveTex:
        tikz_save(Name+'_StrainFiber.tex', show_info=False, strict=False,
                  figureheight='\\figureheight', figurewidth='\\figurewidth',
                  extra_axis_parameters={"axis lines*=left"})
    if SavePng:
        plt.savefig(Name+"_StrainFiber.png", format="png")
    if SaveSvg:
        plt.savefig(Name+"_StrainFiber.svg", format="svg")
    if Show:
        plt.show()


def defineColorMap(ColorMap=0):
    if isinstance(ColorMap, int):
        mapName = ["tab10", "Accent", "Set1_r", "Pastel1", "Set2", "Dark2",
                   "Set3", "jet"]
        nColors = [10, 8, 9, 9, 8, 8, 12, 7]
        cmap = cm.get_cmap(mapName[ColorMap])
        colorList = cmap(np.linspace(0, 1, nColors[ColorMap]))
    elif ColorMap == "unibz":
        colorList = [(0/255, 123/255, 228/255),  (249/255, 107/255, 7/255),
                     (26/255, 140/255, 100/255),  (237/255, 51/255, 26/255),
                     (108/255, 195/255, 171/255), (172/255, 20/255, 90/255)]
    elif ColorMap == "TUM":
        colorList = [(0/255, 101/255, 189/255), (162/255, 173/255, 0/255),
                     (227/255, 114/255, 34/255), (152/255, 198/255, 234/255),
                     (0/255, 82/255, 147/255)]
    elif ColorMap == " A":
        colorList = ["#14325C", "#5398D9", "#F4E3B1", "#D96B0C", "#A53A3B"]
    elif ColorMap == "B":
        colorList = ["#192E5B", "#1D65A6", "#72A2C0", "#00743F", "#F2A104"]
    elif ColorMap == "D":
        colorList = ["#ef476f", "#ffd166", "#06d6a0", "#118ab2", "#073b4c"]
    return(colorList)


if __name__ == "__main__":
    print("test example")
    Laminate1 = laminate()
    Laminate1.SymmetricLaminate()
    Laminate1.calcRuleOfMixtures()
    Laminate1.deg2rad()
    Laminate1.calcPlyPositions()
    Laminate1.calcPoisson21()
    Laminate1.calcABD()
    Laminate1.cleanABD()
    Laminate1.calcInverseABD()
    Laminate1.calcLaminateParameters()
    Laminate1.calcStrain()
    Laminate1.calcStress()
    Laminate1.calcMisesStress()
    Laminate1.calcFailureTsaiWu()
    Laminate1.calcMass()
    print('Laminate mass for one m^2 [t]')
    print(Laminate1.LaminateMass)
    print()
    print("reserve factor after Tsai-Wu:")
    for i in range(Laminate1.nPly):
        print(Laminate1.ReserveTsaiWu[i])
    print("failure mode after Tsai-Wu:")
    for i in range(Laminate1.nPly):
        print(Laminate1.FailureModeTsaiWu[i])
    Laminate1.plotStackLayup()
    Laminate1.plotPlyStrainStress()
    Laminate1.plotPlyStrainStressPly()
    angles = np.rad2deg(Laminate1.theta)
    TW_Res = Laminate1.ReserveTsaiWu
    VMS = Laminate1.stressMisesTop
    Ep = np.zeros(Laminate1.nPly)
    strain = np.zeros(Laminate1.nPly)
    for ii in range(Laminate1.nPly):
        Ep[ii] = np.array([Laminate1.strainLTopVec[ii][0]])
        strain[ii] = np.array([Laminate1.strainPlyTopVec[ii][0]])
#    plotParameterStudyAngle(angles, TW_Res, VMS, Ep, strain)
    Laminate1.writeTablePly(LaminateCooSys=True, PlyCooSys=True,
                            PlyTop=True, PlyBottom=True)
