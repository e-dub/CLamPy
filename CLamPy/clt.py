import numpy as np
import math as mt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Rectangle, Polygon, Patch
import pandas as pd
from tikzplotlib import save as tikz_save


# Ignore NumPy errors for divsion by zero
np.seterr(divide='ignore', invalid='ignore')


class laminate(object):
    """
    A class used to represent a laminate of of fiber-reinforced material.

    Methods
    -------


    Attributes
    ----------
    LaminateName : str
        Name of laminate
    Plot : bool
        Turn on plotting functionality
    nPly : int
        Number of plys in laminate (half the number if symmetrical laminate),
        each value as float
    symmetric : bool
        Declares that laminate is of symmetric plies
    thkPly : list
        Thickness of each ply, each value as float
    theta = list
        Fiber angle of each play, each value as float
    E1 : float
        Young's (elasticity) modulus of ply in fiber direction
    E2 : float
        Young's (elasticity) modulus of ply orthogonal to fiber direction
    G12 : float
        Shear modulus of ply
    nu12 : float
        Poisson's (transverse contration) ratio of ply
    phi : float
        Fiber concentration in percent
    EFiber : float
        Young's (elasticity) modulus of fiber
    EMatrix : float
        Young's (elasticity) modulus of matrix
    nu12Fiber : float
        Poisson's (transverse contration) ratio of fiber
    nu12Matrix : float
        Poisson's (transverse contration) ratio of matrix
    rhoFiber : float
        Density of fiber
    rhoMatrix : float
        Density of fiber
    Nx : float
        Force in laminate x-axis
    Ny : float
        Force in laminate y-axis
    Nxy : float
        Shear force in laminate xy-axis
    Mx : float
        Moment abou laminate x-axis
    My : float
        Moment about laminate y-axis
    Mxy : float
        Moment about laminate xy-axis
    T1 : float
        Tensile strength in fiber direction
    T2 : float
        Tensile strength orthogonal to fiber direction
    C1 : float
        Compressive strength in fiber direction
    C2 : float
        Compressive strength orthongoal to fiber direction
    S12  : float
        Shear strength in ply 12-axis
    S13 : float
        Shear strength in ply 13-axis
    S23 : float
        Shear strength in ply 23-axis
    unitAngles : str {"deg", "rad"}
        Define what unit used for angle measurements

    Todos
    TODO 1: Check and unify units. Should be so done that any unit system of
            consitent units is accepted.  This is not the case.  See mass
            calculation
    """

    LaminateName = ''
    Plot = True
    nPly = 3
    thkPly = [1, 1, 1]
    theta = [-45, -30, 0]
    E1 = 141000
    E2 = 8500
    G12 = 3750
    nu12 = 0.25
    phi = 0.6
    EFiber = 230000
    EMatrix = 3500
    nu12Fiber = 0.2
    nu12Matrix = 0.33
    symmetric = True
    rhoFiber = 1.8
    rhoMatrix = 1.2
    Nx = 100
    Ny = 0
    Nxy = 0
    Mx = 0
    My = 0
    Mxy = 0
    T1 = 1950
    T2 = 48
    C1 = 1480
    C2 = 200
    S12 = 55
    S13 = 55
    S23 = 55
    unitAngles = 'deg'

    def calcRuleOfMixtures(self):
        """
        Calculate ply properties based on fiber and matrix properties.

        Parameters
        ----------
        EFiber : float
            Young's (elasticity) modulus of fiber
        EMatrix : float
            Young's (elasticity) modulus of matrix (assumed to be isotropic)
        E2 : float
            Young's (elasticity) modulus of ply orthogonal to fiber direction
        E1 : float
            Young's (elasticity) modulus of ply in fiber direction
        E2 : float
            Young's (elasticity) modulus of ply orthogonal to fiber direction
        G12Fiber : float
            Shear modulus of fiber
        G12Matrix : float
            Shear modulus of matrix
        G12 : float
            Shear modulus of ply
        nu12 : float
            Poisson's (transverse contration) ratio of ply
        """
        self.E1 = self.EFiber * self.phi + self.EMatrix * (1 - self.phi)
        self.E2 = (
            self.phi / self.EFiber + (1 - self.phi) / self.EMatrix
        ) ** -1
        self.G12Matrix = self.EMatrix / (2 * (1 + self.nu12Matrix))
        self.G12Fiber = self.EFiber / (2 * (1 + self.nu12Fiber))
        self.G12 = (
            self.phi / self.G12Fiber + (1 - self.phi) / self.G12Matrix
        ) ** -1
        # self.G12 = self.G12Matrix*self.G12Fiber/((self.G12Matrix*self.phi)+
        #                                         (1-self.phi)*self.G12Fiber)
        self.nu12 = (
            self.nu12Fiber * self.phi + (1 - self.phi) * self.nu12Matrix
        )

    def SymmetricLaminate(self):
        """
        Mirror plies for symmetric laminates.
        """
        if self.symmetric:
            self.thkPly = np.concatenate(
                (np.array(self.thkPly), np.flip(np.array(self.thkPly))),
                axis=None,
            )
            self.theta = np.concatenate(
                (np.array(self.theta), np.flip(-np.array(self.theta))),
                axis=None,
            )
            self.nPly *= 2

    def calcPlyPositions(self):
        """
        Calculate ply heights and positions in laminate z-axis.

        Parameters
        ----------
        thkPly : float
            Thickness of laminate
        h : list
            Height of each ply, each value is a float
        z : list
            Position of each ply, each value is float
        """
        self.thkLam = sum(self.thkPly)
        self.z = [[]] * self.nPly
        self.h = [[]] * self.nPly
        for i in range(self.nPly):
            self.z[i] = np.array(
                [
                    self.thkLam / 2 - sum(self.thkPly[0:i]),
                    -self.thkLam / 2 + sum(self.thkPly[i + 1 :]),
                ]
            )
            self.h[i] = self.z[i][
                np.argmax((np.abs([self.z[i][0], self.z[i][1]])))
            ]

    def calcPoisson21(self):
        """
        Calculate the Poisson ratio in ply 21-axis.
        """
        self.nu21 = (self.nu12 * self.E2) / self.E1

    def deg2rad(self):
        """
        Set angle measures to radians.
        """
        if self.unitAngles == 'deg':
            self.theta = np.deg2rad(self.theta)
            # self.unitAngles = "rad"

    def calcABD(self):
        """
        Calculate ABD matrix of laminate.
        """
        Q11 = self.E1 / (1 - (self.E2 * self.nu12 ** 2) / self.E1)
        Q22 = self.E2 / self.E1 * Q11
        Q12 = self.nu12 * Q22
        Q66 = self.G12
        Q = np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])
        self.A = np.zeros((3, 3))
        self.B = np.zeros((3, 3))
        self.D = np.zeros((3, 3))
        self.Qbar = [[]] * self.nPly
        self.T = [[]] * self.nPly
        self.Tinv = [[]] * self.nPly
        for i in range(self.nPly):
            c = np.cos((self.theta[i]))
            s = np.sin((self.theta[i]))
            cs = c * s
            cc = c ** 2
            ss = s ** 2
            self.T[i] = np.array(
                [[cc, ss, 2 * cs], [ss, cc, -2 * cs], [-cs, cs, cc - ss]]
            )
            self.Tinv[i] = np.linalg.inv(self.T[i])
            self.Qbar[i] = self.Tinv[i] @ Q @ self.Tinv[i].T
            self.A += self.Qbar[i] * (self.z[i][0] - self.z[i][1])
            self.B += (
                1 / 2 * self.Qbar[i] * (self.z[i][0] ** 2 - self.z[i][1] ** 2)
            )
            self.D += (
                1 / 3 * self.Qbar[i] * (self.z[i][0] ** 3 - self.z[i][1] ** 3)
            )
        self.ABD = np.block([[self.A, self.B], [self.B, self.D]])

    def cleanABD(self):
        """
        Remove numerical errors in the ABD matrix.

        This is carried out by placing small values of the ABD matrix to
        zero.
        """
        self.ABD[abs(self.ABD) < 1e-3] = 0.0

    def calcInverseABD(self):
        """
        Invert ABD matrix.
        """
        self.ABDInv = np.linalg.inv(self.ABD)

    def calcLaminateParameters(self):
        """
        Calculate laminate properties from ABD matrix.
        """
        self.ELx = 1 / (self.ABDInv[0, 0] * self.thkLam)
        self.ELy = 1 / (self.ABDInv[1, 1] * self.thkLam)
        self.GLxy = 1 / (self.ABDInv[2, 2] * self.thkLam)
        self.nuLxy = self.ABDInv[0, 1] / self.ABDInv[0, 0]
        self.nuLyx = self.ABDInv[0, 1] / self.ABDInv[1, 1]

    def calcStrain(self):
        """
        Calculate strains in laminate from load.
        """
        loadVec = np.array(
            [[self.Nx, self.Ny, self.Nxy, self.Mx, self.My, self.Mxy]]
        ).T
        self.strainVec = self.ABDInv @ loadVec
        # self.strainLamVec = [[]]*self.nPly
        self.strainLamTopVec = [[]] * self.nPly
        self.strainLamBotVec = [[]] * self.nPly
        self.strainPlyTopVec = [[]] * self.nPly
        self.strainPlyBotVec = [[]] * self.nPly
        for i in range(self.nPly):
            if self.z[i][0] > 0:
                self.strainLamTopVec[i] = self.strainVec[:3] + self.strainVec[
                    3:
                ] * (self.z[i][0])
                self.strainLamBotVec[i] = self.strainVec[:3] + self.strainVec[
                    3:
                ] * (self.z[i][1])
            else:
                self.strainLamTopVec[i] = self.strainVec[:3] + self.strainVec[
                    3:
                ] * (self.z[i][1])
                self.strainLamBotVec[i] = self.strainVec[:3] + self.strainVec[
                    3:
                ] * (self.z[i][0])
            self.strainPlyTopVec[i] = (
                np.linalg.inv(self.T[i].T) @ self.strainLamTopVec[i]
            )
            self.strainPlyBotVec[i] = (
                np.linalg.inv(self.T[i].T) @ self.strainLamBotVec[i]
            )

    def calcStress(self):
        """
        Calculate stress components from strains.
        """
        self.stressLamTopVec = [[]] * self.nPly
        self.stressLamBotVec = [[]] * self.nPly
        self.stressPlyTopVec = [[]] * self.nPly
        self.stressPlyBotVec = [[]] * self.nPly
        for i in range(self.nPly):
            self.stressLamTopVec[i] = self.Qbar[i] @ self.strainLamTopVec[i]
            self.stressLamBotVec[i] = self.Qbar[i] @ self.strainLamBotVec[i]
            self.stressPlyTopVec[i] = self.T[i] @ self.stressLamTopVec[i]
            self.stressPlyBotVec[i] = self.T[i] @ self.stressLamTopVec[i]

    def calcMisesStress(self):
        """
        Calculate equivalent stress after von Mises from stress components.
        """
        A = np.array([[1.0, -0.5, 0.0], [-0.5, 1.0, 0.0], [0.0, 0.0, 3.0]])
        self.stressMisesTop = np.zeros((self.nPly,))
        self.stressMisesBot = np.zeros((self.nPly,))
        for i in range(self.nPly):
            self.stressMisesBot[i] = np.sqrt(
                self.stressLamBotVec[i].T @ A @ self.stressLamBotVec[i]
            )
            self.stressMisesTop[i] = np.sqrt(
                self.stressLamTopVec[i].T @ A @ self.stressLamTopVec[i]
            )

    def calcFailureTsaiWu(self):
        """
        Calculate Tsai-Wu failure criterion from stress components.
        """
        self.FailureTsaiWu = []
        self.ReserveTsaiWu = []
        F11 = 1 / (self.T1 * self.C1)
        F22 = 1 / (self.T2 * self.C2)
        F1 = (1 / self.T1) - (1 / self.C1)
        F2 = (1 / self.T2) - (1 / self.C2)
        F12 = -0.5 / np.sqrt(self.T1 * self.T2 * self.C1 * self.C2)
        F33 = 1 / (self.S12 ** 2)
        self.FailureModeTsaiWu = [[]] * self.nPly
        self.FailureTsaiWu = [[]] * self.nPly
        self.ReserveTsaiWu = [[]] * self.nPly
        for i in range(self.nPly):
            Mode = [[]] * 3
            aTop = (
                F11 * self.stressPlyTopVec[i][0] ** 2
                + F22 * self.stressPlyTopVec[i][1] ** 2
                + F33 * self.stressPlyTopVec[i][2] ** 2
                + 2
                * F12
                * self.stressPlyTopVec[i][0]
                * self.stressPlyTopVec[i][1]
            )[0]
            aBot = (
                F11 * self.stressPlyBotVec[i][0] ** 2
                + F22 * self.stressPlyBotVec[i][1] ** 2
                + F33 * self.stressPlyBotVec[i][2] ** 2
                + 2
                * F12
                * self.stressPlyBotVec[i][0]
                * self.stressPlyBotVec[i][1]
            )[0]
            bTop = (
                F1 * self.stressPlyTopVec[i][0]
                + F2 * self.stressPlyTopVec[i][1]
            )[0]
            bBot = (
                F1 * self.stressPlyBotVec[i][0]
                + F2 * self.stressPlyBotVec[i][1]
            )[0]
            a = max(aTop, aBot)
            b = max(bTop, bBot)
            self.ReserveTsaiWu[i] = (-b + (b ** 2 + 4 * a) ** (1 / 2)) / (
                2 * a
            )
            self.FailureTsaiWu[i] = (2 * a) / (
                -b + (b ** 2 + 4 * a) ** (1 / 2)
            )
            if self.ReserveTsaiWu[i] < 1.0:
                Modes = ['fiber', 'matrix', 'shear']
                Mode[0] = max(
                    abs(
                        F1 * self.stressPlyTopVec[i][0]
                        + F11 * self.stressPlyTopVec[i][0] ** 2
                    ),
                    abs(
                        F1 * self.stressPlyBotVec[i][0]
                        + F11 * self.stressPlyBotVec[i][0] ** 2
                    ),
                )
                Mode[1] = max(
                    abs(
                        F2 * self.stressPlyTopVec[i][1]
                        + F22 * self.stressPlyTopVec[i][1] ** 2
                    ),
                    abs(
                        F2 * self.stressPlyBotVec[i][1]
                        + F22 * self.stressPlyBotVec[i][1] ** 2
                    ),
                )
                Mode[2] = max(
                    abs(F33 * self.stressPlyTopVec[i][2] ** 2),
                    abs(F33 * self.stressPlyBotVec[i][2] ** 2),
                )
                self.FailureModeTsaiWu[i] = Modes[np.argmax(Mode)]
            else:
                self.FailureModeTsaiWu[i] = 'none'

    def calcMass(self):
        """
        Calculate laminate mass.
        """
        self.LaminateVolume = (
            np.sum(self.thkPly) / 1000
        )  # fibers volume in m^3 for one m^2 of laminate
        self.LaminateMass = self.LaminateVolume * (
            self.rhoFiber * self.phi + self.rhoMatrix * (1 - self.phi)
        )
        # Laminate mass in tons, referred to one m^2

    def plotStackLayup(
        self,
        savePng=True,
        saveSvg=True,
        saveTex=True,
        plotShow=True,
        fontSize=12,
        fontName='Tex Gyre Pagella',
        ColorMap=0,
    ):
        """
        Plot ply stackup of laminate.
        """
        # if fontName == "Palatino":
        #    fontName = "Tex Gyre Pagella"
        colorList = defineColorMap(
            ColorMap, self.nPly - int(self.symmetric * self.nPly / 2)
        )
        fig = plt.figure(
            figsize=(2 + int(self.nPly - 1) / 2, 2 + int(self.nPly - 1) / 2)
        )
        ax = fig.add_subplot()

        def plot_hatches(ax, angle, ii, delta, offset=0.0):
            nPoints = 501
            xL = -2
            xU = 2
            x = np.linspace(xL, xU, nPoints)
            for c in np.arange(xL, xU, offset):
                yprime = mt.cos(angle) * c + mt.sin(angle) * x + delta
                xprime = -mt.sin(angle) * c + mt.cos(angle) * x + delta
                yprime[yprime > (0.999 + delta)] = np.nan
                xprime[xprime > (0.999 + delta)] = np.nan
                yprime[yprime < (0.001 + delta)] = np.nan
                yprime[yprime < (0.001 + delta)] = np.nan
                xprime[xprime < (0.001 + delta)] = np.nan
                xprimeTemp = xprime[np.logical_not(np.isnan(xprime))]
                yprimeTemp = yprime[np.logical_not(np.isnan(xprime))]
                xprime = xprimeTemp[np.logical_not(np.isnan(yprimeTemp))]
                yprime = yprimeTemp[np.logical_not(np.isnan(yprimeTemp))]
                ax.plot(xprime, yprime, color='k', zorder=ii)
            return ax

        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly - 1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = i
            index[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = -1
        colorIndex -= min(colorIndex)
        for ii in range(self.nPly):
            el = ii * 0.25
            plot_hatches(ax, angle[ii], ii, el, offset=0.05)
            ax.add_patch(
                Rectangle(
                    (0 + el, 0 + el),
                    1,
                    1,
                    fill='r',
                    alpha=0.85,
                    color=colorList[int(colorIndex[ii])],
                    zorder=ii,
                )
            )
            plt.text(
                1.05 + el,
                0.05 + el,
                str(round(formatNumber(round(angle[ii], 4)), 2)) + ' deg',
                fontsize=fontSize,
            )
        ax.set_ylim(0, 1 + el)
        ax.set_xlim(0, 1 + el)
        plt.axis('off')
        # plt.rcParams["font.family"] = fontName
        if self.LaminateName == '':
            Name = ''
        else:
            Name = self.LaminateName + '_'
        if savePng:
            plt.savefig(
                Name + 'StackLayup.png', format='png', bbox_inches='tight'
            )
        if saveSvg:
            plt.savefig(
                Name + 'StackLayup.svg', format='svg', bbox_inches='tight'
            )
        if saveTex:  # Currently not working correctly, patches not shown.
            tikz_save(
                Name + 'StackLayup.tex',
                show_info=False,
                strict=False,
                extra_axis_parameters={'axis lines*=left'},
            )
        if plotShow:
            plt.show()

    def plotPlyStrainStress(
        self,
        Show=True,
        SaveTex=True,
        SavePng=True,
        SaveSvg=True,
        ColorMap=0,
        Legend=True,
    ):
        """
        Plot stresses and strains of laminate.
        """
        colorList = defineColorMap(
            ColorMap, self.nPly - int(self.symmetric * self.nPly / 2)
        )
        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly - 1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = i
            index[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = -1
        colorIndex -= min(colorIndex)
        legendElements = []
        angleList = []
        angleOrder = np.argsort(angle)
        for i in range(self.nPly):
            j = angleOrder[i]
            val = abs(self.theta[j])
            if val not in angleList:
                angleList.append(val)
                lab = (
                    str(round(formatNumber(round(abs(angle[j]), 4)), 2))
                    + ' deg'
                )
                legendElements.append(
                    Patch(
                        facecolor=colorList[int(colorIndex[j])],
                        edgecolor=colorList[int(colorIndex[j])],
                        label=lab,
                    )
                )
        fig = plt.figure(figsize=(7, 3))
        ax = fig.add_subplot()
        width = 1.5
        xMin = -width / 2
        Shift = 2
        xShift = []
        x = [[]] * self.nPly
        for i in range(self.nPly):
            ax.add_patch(
                Rectangle(
                    (xMin, self.z[i][0]),
                    width,
                    -self.thkPly[i],
                    fill='r',
                    alpha=0.85,
                    color=colorList[int(colorIndex[i])],
                    zorder=i,
                )
            )
        for ii in range(3):
            xShift.append(1 * Shift + 2 * Shift * ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (
                        np.array(
                            [
                                min(self.strainLamTopVec[i][ii][0], 0),
                                max(self.strainLamTopVec[i][ii][0], 0),
                                max(self.strainLamBotVec[i][ii][0], 0),
                                min(self.strainLamBotVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.strainLamTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
                else:
                    x[i] = (
                        np.array(
                            [
                                max(self.strainLamBotVec[i][ii][0], 0),
                                min(self.strainLamBotVec[i][ii][0], 0),
                                min(self.strainLamTopVec[i][ii][0], 0),
                                max(self.strainLamTopVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.strainLamTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
            xShift.append(2 * Shift + 2 * Shift * ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (
                        np.array(
                            [
                                min(self.stressLamTopVec[i][ii][0], 0),
                                max(self.stressLamTopVec[i][ii][0], 0),
                                max(self.stressLamBotVec[i][ii][0], 0),
                                min(self.stressLamBotVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.stressLamTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
                else:
                    x[i] = (
                        np.array(
                            [
                                min(self.stressLamBotVec[i][ii][0], 0),
                                max(self.stressLamBotVec[i][ii][0], 0),
                                max(self.stressLamTopVec[i][ii][0], 0),
                                min(self.stressLamTopVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.stressLamTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
        plt.ylabel('laminate position $z$ [mm]')
        ax.set_ylim([np.min(self.z), np.max(self.z)])
        # plt.xticks(xShift, ('0\n strain x', '0\n stress x', '0\n strain y', '0\n stress y',
        #                    '0\nstrain xy', '0\n stress xy'))
        # plt.xticks(xShift, ['0', '0', '0','0', '0', '0'])
        # HERE AN ATTEMPT TO MAKE LABELS WORK IN TIKZ...
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_label_position('bottom')
        ax.set_xticks(xShift)
        ax.set_xticklabels(('0', '0', '0', '0', '0', '0'))
        plt.xlabel('values normalized based on maximum [-]')
        ax2 = ax.twiny()
        ax2.xaxis.set_tick_params(labeltop='on')
        # ax2.xaxis.set_ticks_position("top")
        # ax2.xaxis.set_label_position("top")
        # ax2.spines["top"].set_position(("axes", -0.25))
        # ax2.set_frame_on(False)
        # ax2.patch.set_visible(False)
        # ax2.spines["top"].set_visible(False)
        ax2.set_xticks(np.insert(xShift, 0, 0))
        ax2.set_xticklabels(
            (
                'stackup',
                'strain x',
                'stress x',
                'strain y',
                'stress y',
                'strain xy',
                'stress xy',
            )
        )
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
            plt.legend(
                handles=legendElements[::-1],
                loc='center left',
                bbox_to_anchor=(1.05, 0.5),
                frameon=False,
            )
        if SaveTex:
            tikz_save(
                'LaminateStrainStress.tex', show_info=False, strict=False
            )  # ,figureheight='\\figureheight', figurewidth='\\figurewidth', extra_axis_parameters={"axis lines*=left"})
        if SavePng:
            plt.savefig('LaminateStrainStress.png', format='png')
        if SaveSvg:
            plt.savefig('LaminateStrainStress.svg', format='svg')
        if Show:
            plt.show()

    def plotPlyStrainStressPly(
        self,
        Show=True,
        SaveTex=True,
        SavePng=True,
        SaveSvg=True,
        ColorMap=0,
        Legend=True,
    ):
        """
        Plot strain and stress in ply.
        """
        colorList = defineColorMap(
            ColorMap, self.nPly - int(self.symmetric * self.nPly / 2)
        )
        angle = np.rad2deg(self.theta)
        index = np.linspace(0, self.nPly - 1, self.nPly)
        colorIndex = np.ones((self.nPly))
        for i in range(self.nPly):
            if max(index) < 0:
                break
            colorIndex[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = i
            index[
                (np.abs(np.abs(angle) - np.abs(angle[i])) <= 0.001).tolist()
            ] = -1
        colorIndex -= min(colorIndex)
        legendElements = []
        angleList = []
        angleOrder = np.argsort(angle)
        for i in range(self.nPly):
            j = angleOrder[i]
            val = abs(self.theta[j])
            if val not in angleList:
                angleList.append(val)
                lab = (
                    str(round(formatNumber(round(abs(angle[j]), 4)), 2))
                    + ' deg'
                )
                legendElements.append(
                    Patch(
                        facecolor=colorList[int(colorIndex[j])],
                        edgecolor=colorList[int(colorIndex[j])],
                        label=lab,
                    )
                )
        fig = plt.figure(figsize=(7, 3))
        ax = fig.add_subplot()
        width = 1.5
        xMin = -width / 2
        Shift = 2
        xShift = []
        x = [[]] * self.nPly
        for i in range(self.nPly):
            ax.add_patch(
                Rectangle(
                    (xMin, self.z[i][0]),
                    width,
                    -self.thkPly[i],
                    fill='r',
                    alpha=0.85,
                    color=colorList[int(colorIndex[i])],
                    zorder=i,
                )
            )
        for ii in range(3):
            xShift.append(1 * Shift + 2 * Shift * ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (
                        np.array(
                            [
                                min(self.strainPlyTopVec[i][ii][0], 0),
                                max(self.strainPlyTopVec[i][ii][0], 0),
                                max(self.strainPlyBotVec[i][ii][0], 0),
                                min(self.strainPlyBotVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.strainPlyTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
                else:
                    x[i] = (
                        np.array(
                            [
                                max(self.strainPlyBotVec[i][ii][0], 0),
                                min(self.strainPlyBotVec[i][ii][0], 0),
                                min(self.strainPlyTopVec[i][ii][0], 0),
                                max(self.strainPlyTopVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.strainPlyTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
            xShift.append(2 * Shift + 2 * Shift * ii)
            for i in range(self.nPly):
                if self.z[i][0] > 0:
                    x[i] = (
                        np.array(
                            [
                                min(self.stressPlyTopVec[i][ii][0], 0),
                                max(self.stressPlyTopVec[i][ii][0], 0),
                                max(self.stressPlyBotVec[i][ii][0], 0),
                                min(self.stressPlyBotVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.stressPlyTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
                else:
                    x[i] = (
                        np.array(
                            [
                                min(self.stressPlyBotVec[i][ii][0], 0),
                                max(self.stressPlyBotVec[i][ii][0], 0),
                                max(self.stressPlyTopVec[i][ii][0], 0),
                                min(self.stressPlyTopVec[i][ii][0], 0),
                            ]
                        )
                        / np.max(np.max(self.stressPlyTopVec))
                        + xShift[-1]
                    )
                    y = np.array(
                        [
                            self.z[i][0],
                            self.z[i][0],
                            self.z[i][1],
                            self.z[i][1],
                        ]
                    )
                    ax.add_patch(
                        Polygon(
                            xy=list(zip(x[i], y)),
                            fill='r',
                            alpha=0.85,
                            color=colorList[int(colorIndex[i])],
                            zorder=i,
                        )
                    )
        plt.ylabel('laminate position $z$ [mm]')
        ax.set_ylim([np.min(self.z), np.max(self.z)])
        # plt.xticks(xShift, ('0\n strain x', '0\n stress x', '0\n strain y', '0\n stress y',
        #                    '0\nstrain xy', '0\n stress xy'))
        # plt.xticks(xShift, ['0', '0', '0','0', '0', '0'])
        # HERE AN ATTEMPT TO MAKE LABELS WORK IN TIKZ...
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_label_position('bottom')
        ax.set_xticks(xShift)
        ax.set_xticklabels(('0', '0', '0', '0', '0', '0'))
        plt.xlabel('values normalized based on maximum [-]')
        ax2 = ax.twiny()
        ax2.xaxis.set_tick_params(labeltop='on')
        # ax2.xaxis.set_ticks_position("top")
        # ax2.xaxis.set_label_position("top")
        # ax2.spines["top"].set_position(("axes", -0.25))
        # ax2.set_frame_on(False)
        # ax2.patch.set_visible(False)
        # ax2.spines["top"].set_visible(False)
        ax2.set_xticks(np.insert(xShift, 0, 0))
        ax2.set_xticklabels(
            (
                'stackup',
                'strain 1',
                'stress 1',
                'strain 2',
                'stress 2',
                'strain 12',
                'stress 12',
            )
        )
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
            plt.legend(
                handles=legendElements[::-1],
                loc='center left',
                bbox_to_anchor=(1.05, 0.5),
                frameon=False,
            )
        if SaveTex:
            tikz_save('PlyStrainStress.tex', show_info=False, strict=False)
            # figureheight='\\figureheight', figurewidth='\\figurewidth',
            # extra_axis_parameters={"axis lines*=left"})
        if SavePng:
            plt.savefig('PlyStrainStress.png', format='png')
        if SaveSvg:
            plt.savefig('PlyStrainStress.svg', format='svg')
        if Show:
            plt.show()

    def writeTablePly(
        self,
        writeTex=True,
        writeCsv=True,
        printData=True,
        DecimalPlaces=4,
        DecimalSymbol='.',
        LaminateCS=True,
        PlyCS=True,
        PlyTop=True,
        PlyBottom=True,
    ):
        """
        Write table of ply values.
        """
        FileName = 'PlyTable'
        TestDataTop = {
            'ply': np.arange(1, self.nPly + 1),
            'fiber orientation': np.rad2deg(self.theta),
            'strain x': [i[0][0] for i in self.strainLamTopVec],
            'stress x': [i[0][0] for i in self.stressLamTopVec],
            'strain y': [i[1][0] for i in self.strainLamTopVec],
            'stress y': [i[1][0] for i in self.stressLamTopVec],
            'strain xy': [i[2][0] for i in self.strainLamTopVec],
            'stress xy': [i[2][0] for i in self.stressLamTopVec],
            'strain 1': [i[0][0] for i in self.strainPlyTopVec],
            'stress 1': [i[0][0] for i in self.stressPlyTopVec],
            'strain 2': [i[1][0] for i in self.strainPlyTopVec],
            'stress 2': [i[1][0] for i in self.stressPlyTopVec],
            'strain 12': [i[2][0] for i in self.strainPlyTopVec],
            'stress 12': [i[2][0] for i in self.stressPlyTopVec],
            'Tsai-Wu failure': self.FailureTsaiWu,
            'Tsai-Wu reserve': self.ReserveTsaiWu,
        }
        TestDataBot = {
            'ply': np.arange(1, self.nPly + 1),
            'fiber orientation': np.rad2deg(self.theta),
            'strain x': [i[0][0] for i in self.strainLamBotVec],
            'stress x': [i[0][0] for i in self.stressLamBotVec],
            'strain y': [i[1][0] for i in self.strainLamBotVec],
            'stress y': [i[1][0] for i in self.stressLamBotVec],
            'strain xy': [i[2][0] for i in self.strainLamBotVec],
            'stress xy': [i[2][0] for i in self.stressLamBotVec],
            'strain 1': [i[0][0] for i in self.strainPlyBotVec],
            'stress 1': [i[0][0] for i in self.stressPlyBotVec],
            'strain 2': [i[1][0] for i in self.strainPlyBotVec],
            'stress 2': [i[1][0] for i in self.stressPlyBotVec],
            'strain 12': [i[2][0] for i in self.strainPlyBotVec],
            'stress 12': [i[2][0] for i in self.stressPlyBotVec],
            'Tsai-Wu failure': self.FailureTsaiWu,
            'Tsai-Wu reserve': self.ReserveTsaiWu,
        }
        ColumnLam = [
            'ply',
            'fiber orientation',
            'strain x',
            'stress x',
            'strain y',
            'stress y',
            'strain xy',
            'stress xy',
            'Tsai-Wu failure',
            'Tsai-Wu reserve',
        ]
        ColumnPly = [
            'ply',
            'fiber orientation',
            'strain 1',
            'stress 1',
            'strain 2',
            'stress 2',
            'strain 12',
            'stress 12',
            'Tsai-Wu failure',
            'Tsai-Wu reserve',
        ]
        LaminateDataFrameTop = pd.DataFrame(TestDataTop, columns=ColumnLam)
        LaminateDataFrameBot = pd.DataFrame(TestDataBot, columns=ColumnLam)
        PlyDataFrameTop = pd.DataFrame(TestDataTop, columns=ColumnPly)
        PlyDataFrameBot = pd.DataFrame(TestDataBot, columns=ColumnPly)
        LaminateDataFrameTop = LaminateDataFrameTop.round(DecimalPlaces)
        LaminateDataFrameBot = LaminateDataFrameBot.round(DecimalPlaces)
        PlyDataFrameTop = PlyDataFrameTop.round(DecimalPlaces)
        PlyDataFrameBot = PlyDataFrameBot.round(DecimalPlaces)
        if writeTex:
            if LaminateCS:
                headerTex = [
                    'ply',
                    'fiber orientation $\\theta$ [deg]',
                    'strain x $\\varepsilon_x$ [-]',
                    'stress x $\\sigma_x$ [MPa]',
                    'strain y $\\varepsilon_y$ [-]',
                    'stress y $\\sigma_y$ [MPa]',
                    'strain xy $\\varepsilon_{xy}$ [-]',
                    'stress xy $\\sigma_{xy}$ [MPa]',
                    'Tsai-Wu failure $F_\\mathrm{TW}$ [-]',
                    'Tsai-Wu reserve $R_\\mathrm{TW}$ [-]',
                ]
                if PlyTop:
                    LaminateDataFrameTop.to_latex(
                        FileName + 'LaminateTop.tex',
                        index=None,
                        decimal=DecimalSymbol,
                        header=headerTex,
                        escape=False,
                    )
                    if writeCsv:
                        LaminateDataFrameTop.to_excel(
                            FileName + 'LaminateTop.xlsx',
                            index=None,
                            header=True,
                        )
                    if printData:
                        print(
                            self.LaminateName,
                            'laminate properties in top reference system:',
                        )
                        print(LaminateDataFrameTop.to_string(index=False))

                if PlyBottom:
                    LaminateDataFrameBot.to_latex(
                        FileName + 'LaminateBot.tex',
                        index=None,
                        decimal=DecimalSymbol,
                        header=headerTex,
                        escape=False,
                    )
                    if writeCsv:
                        LaminateDataFrameBot.to_excel(
                            FileName + 'LaminateBot.xlsx',
                            index=None,
                            header=True,
                        )
                    if printData:
                        print(
                            self.LaminateName,
                            'laminate properties in bottom reference system:',
                        )
                        print(LaminateDataFrameBot.to_string(index=False))
            if PlyCS:
                headerTex = [
                    'ply',
                    'fiber orientation $\\theta$ [deg]',
                    'strain 1 $\\varepsilon_1$ [-]',
                    'stress 1 $\\sigma_1$ [MPa]',
                    'strain 2 $\\varepsilon_2$ [-]',
                    'stress 2 $\\sigma_2$ [MPa]',
                    'strain 12 $\\varepsilon_{12}$ [-]',
                    'stress 12 $\\sigma_{12}$ [MPa]',
                    'Tsai-Wu failure $F_\\mathrm{TW}$ [-]',
                    'Tsai-Wu reserve $R_\\mathrm{TW}$ [-]',
                ]
                if PlyTop:
                    PlyDataFrameTop.to_latex(
                        FileName + 'PlyTop.tex',
                        index=None,
                        decimal=DecimalSymbol,
                        header=headerTex,
                        escape=False,
                    )
                    if writeCsv:
                        PlyDataFrameTop.to_excel(
                            FileName + 'PlyTop.xlsx', index=None, header=True
                        )
                    if printData:
                        print(
                            self.LaminateName,
                            'ply properties in top reference system:',
                        )
                        print(PlyDataFrameTop.to_string(index=False))
                if PlyBottom:
                    PlyDataFrameBot.to_latex(
                        FileName + 'PlyBot.tex',
                        index=None,
                        decimal=DecimalSymbol,
                        header=headerTex,
                        escape=False,
                    )
                    if writeCsv:
                        PlyDataFrameBot.to_excel(
                            FileName + 'PlyBot.xlsx', index=None, header=True
                        )
                    if printData:
                        print(
                            self.LaminateName,
                            'ply properties in bottom reference system:',
                        )
                        print(PlyDataFrameBot.to_string(index=False))


def formatNumber(value):
    """
    Format numbers so that floats with no decimals are integers.
    """
    if value % 1 == 0:
        value = int(value)
    return value


def plotParameterStudyAngle(
    angles,
    RTW,
    sigmaM,
    epsilon1,
    epsilonx,
    Show=True,
    SaveTex=True,
    SavePng=True,
    SaveSvg=True,
    Name='',
    buffer=0.1,
    Grid=False,
    plotSize=(7, 5),
    ColorTW='r',
    SymbolTW='.',
    LineTW='-',
    ColorM='b',
    SymbolM='.',
    LineM='-',
    ColorStrainFiber='g',
    SymbolStrainFiber='.',
    LineStrainFiber='-',
    ColorStrainx='g',
    SymbolStrainx='.',
    LineStrainx='-',
    xBuffer=0,
):
    """
    Plot parameter study of ply angles.

    Parameters
    ----------
    angles : list or NumPy array
        Angle values of parameter study
    RTW : list or NumPy array
        Tsai-Wu reserve factor of parameter study
    sigmaM : list or NumPy array
        Mises stress of parameter study
    epsilon1 : list or NumPy array
        ply strain in fiber direction
    epsilonx : list or NumPy array
        ply strain in laminate direction x
    Show : bool
        Show figures, if off only files will be generated
    SaveTex : bool
        Save figure in TeX format using tikz
    SavePng : bool
        Save figure in PNG format raster graphics
    SaveSvg : bool
        Save figure in SVG format vector graphics
    Name : str
        Name of study to be added to file names
    buffer : float
        Add space to plots for formatting in percent of difference between
        maximum and minimum values
    Gird : bool
        Add grid to plot
    plotSize : tuple of two floats
        Define size of plots, especially important for PNGs
    ColorTW : str
        Color for Tsai-Wu plot
    SymbolTW : str
        Symbol for Tsai-Wu plot
    LineTW : str
        Line type for Tsai-Wu plot
    ColorM : str
        Color for Mises stress plot
    SymbolM : str
        Symbol for Mises stress plot
    LineM : str
        Line type for Mises stress plot
    ColorStrain1 : str
        Color for fiber strain plot
    SymbolStrain1 : str
        Symbol for fiber strain plot
    LineStrain1 : str
        Line type for fiber strain plot
    ColorStrainX : str
        Color for laminate strain plot
    SymbolStrainX : str
        Symbol for laminate strain plot
    LineStrainX : str
        Line type for laminate strain plot
    xBuffer : float
        Add space to x-axis in plots
    """
    fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    plt.plot(angles, RTW, ColorTW + SymbolTW + LineTW)
    plt.ylabel('Tsai-Wu reserve factor $R_{TW}$ [-]')
    plt.xlim([min(angles) - xBuffer, max(angles) + xBuffer])
    plt.ylim(
        [
            min(RTW) - max(np.abs(RTW)) * buffer,
            max(RTW) + max(np.abs(RTW)) * buffer,
        ]
    )
    yticks = np.arange(0, max(RTW) + 10, step=10)
    yticks[0] = 1
    plt.yticks(yticks)
    xticks = np.arange(min(angles), max(angles) + 1, step=15)
    plt.xticks(xticks)
    plt.xlabel('fiber orientation angle $\\theta$ [deg]')
    if SaveTex:
        tikz_save(
            Name + '_TsaiWuReserve.tex',
            show_info=False,
            strict=False,
            extra_axis_parameters={
                'axis lines*=left',
                'height=\\figureheight',
                'width=\\figurewidth',
            },
        )
    if SavePng:
        plt.savefig(Name + '_TsaiWuReserve.png', format='png')
    if SaveSvg:
        plt.savefig(Name + '_TsaiWuReserve.svg', format='svg')
    if Show:
        plt.show()
    fig, ax2 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.yaxis.set_ticks_position('left')
    ax2.xaxis.set_ticks_position('bottom')
    plt.plot(angles, sigmaM, ColorM + SymbolM + LineM)
    plt.ylabel('equivalent stress after Mises $\\sigma_M$ [MPa] ')
    plt.xlim([min(angles) - xBuffer, max(angles) + xBuffer])
    plt.ylim(
        [
            min(sigmaM) - max(np.abs(sigmaM)) * buffer,
            max(sigmaM) + max(np.abs(sigmaM)) * buffer,
        ]
    )
    xticks = np.arange(min(angles), max(angles) + 1, step=15)
    plt.xticks(xticks)
    plt.xlabel('fiber orientation angle $\\theta$ [deg]')
    plt.tight_layout()
    if SaveTex:
        tikz_save(
            Name + '_StressEqv.tex',
            show_info=False,
            strict=False,
            extra_axis_parameters={
                'axis lines*=left',
                'height=\\figureheight',
                'width=\\figurewidth',
            },
        )
    if SavePng:
        plt.savefig(Name + '_StressMises.png', format='png')
    if SaveSvg:
        plt.savefig(Name + '_StressMises.svg', format='svg')
    if Show:
        plt.show()
    fig, ax3 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax3.spines['right'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.yaxis.set_ticks_position('left')
    ax3.xaxis.set_ticks_position('bottom')
    plt.plot(angles, epsilon1, ColorStrainx + SymbolStrainx + LineStrainx)
    plt.ylabel('epsilonx in $x$ $\\varepsilon_x$ [-]')
    plt.xlim([min(angles) - xBuffer, max(angles) + xBuffer])
    plt.ylim(
        [
            min(epsilon1) - max(np.abs(epsilon1)) * buffer,
            max(epsilon1) + max(np.abs(epsilon1)) * buffer,
        ]
    )
    xticks = np.arange(min(angles), max(angles) + 1, step=15)
    plt.xticks(xticks)
    plt.xlabel('fiber orientation angle $\\theta$ [deg]')
    plt.tight_layout()
    if SaveTex:
        tikz_save(
            Name + '_Strain_x.tex',
            show_info=False,
            strict=False,
            extra_axis_parameters={
                'axis lines*=left',
                'height=\\figureheight',
                'width=\\figurewidth',
            },
        )
    if SavePng:
        plt.savefig(Name + '_Strain_x.png', format='png')
    if SaveSvg:
        plt.savefig(Name + '_Strain_x.svg', format='svg')
    if Show:
        plt.show()
    fig, ax4 = plt.subplots(nrows=1, ncols=1, figsize=plotSize)
    plt.grid(Grid)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.yaxis.set_ticks_position('left')
    ax4.xaxis.set_ticks_position('bottom')
    plt.plot(
        angles,
        epsilonx,
        ColorStrainFiber + SymbolStrainFiber + LineStrainFiber,
    )
    plt.ylabel('strain in $\\parallel$ $\\varepsilon_{\\parallel}$ [-]')
    plt.xlim([min(angles) - xBuffer, max(angles) + xBuffer])
    plt.ylim(
        [
            min(epsilonx) - max(np.abs(epsilonx)) * buffer,
            max(epsilonx) + max(np.abs(epsilonx)) * buffer,
        ]
    )
    xticks = np.arange(min(angles), max(angles) + 1, step=15)
    plt.xticks(xticks)
    plt.xlabel('fiber orientation angle $\\theta$ [deg]')
    plt.tight_layout()
    if SaveTex:
        tikz_save(
            Name + '_StrainFiber.tex',
            show_info=False,
            strict=False,
            extra_axis_parameters={
                'axis lines*=left',
                'height=\\figureheight',
                'width=\\figurewidth',
            },
        )
    if SavePng:
        plt.savefig(Name + '_StrainFiber.png', format='png')
    if SaveSvg:
        plt.savefig(Name + '_StrainFiber.svg', format='svg')
    if Show:
        plt.show()


def defineColorMap(ColorMap=0, nColors=3):
    """
    Define color map.

    Parameters
    ----------
    ColorMap : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    colorList : NumPy array
        List of colors

    """
    mapName = [
        'tab10',
        'Accent',
        'Set1_r',
        'Pastel1',
        'Set2',
        'Dark2',
        'Set3',
        'jet',
        'tab20',
    ]
    nColorMap = [10, 8, 9, 9, 8, 8, 12, 7, 20]
    if isinstance(ColorMap, int):
        if nColors > 10:
            cmap = cm.get_cmap('tab20')
            colorList = cmap(np.linspace(0, 1, nColorMap[-1]))
            colorList = np.concatenate((colorList[::2], colorList[1::2]), 0)
        else:
            cmap = cm.get_cmap(mapName[ColorMap])
            colorList = cmap(np.linspace(0, 1, nColorMap[ColorMap]))
    elif ColorMap == 'unibz':
        colorList = [
            (0 / 255, 123 / 255, 228 / 255),
            (249 / 255, 107 / 255, 7 / 255),
            (26 / 255, 140 / 255, 100 / 255),
            (237 / 255, 51 / 255, 26 / 255),
            (108 / 255, 195 / 255, 171 / 255),
            (172 / 255, 20 / 255, 90 / 255),
        ]
    elif ColorMap == 'TUM':
        colorList = [
            (0 / 255, 101 / 255, 189 / 255),
            (162 / 255, 173 / 255, 0 / 255),
            (227 / 255, 114 / 255, 34 / 255),
            (152 / 255, 198 / 255, 234 / 255),
            (0 / 255, 82 / 255, 147 / 255),
        ]
    elif ColorMap == ' A':
        colorList = ['#14325C', '#5398D9', '#F4E3B1', '#D96B0C', '#A53A3B']
    elif ColorMap == 'B':
        colorList = ['#192E5B', '#1D65A6', '#72A2C0', '#00743F', '#F2A104']
    elif ColorMap == 'D':
        colorList = ['#ef476f', '#ffd166', '#06d6a0', '#118ab2', '#073b4c']
    return colorList


if __name__ == '__main__':
    print('test example')
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
    print('Laminate mass for 1000 mm² (1 m²) [t]')
    print(Laminate1.LaminateMass)
    print()
    print('reserve factor after Tsai-Wu:')
    for i in range(Laminate1.nPly):
        print(
            '  Ply '
            + str(i)
            + ': '
            + str(round(Laminate1.ReserveTsaiWu[i], 4))
        )
    print('failure mode after Tsai-Wu:')
    for i in range(Laminate1.nPly):
        print('  Ply ' + str(i) + ': ' + Laminate1.FailureModeTsaiWu[i])
    Laminate1.plotStackLayup()
    Laminate1.plotPlyStrainStress()
    Laminate1.plotPlyStrainStressPly()
    angles = np.rad2deg(Laminate1.theta)
    RTW = Laminate1.ReserveTsaiWu
    sigmaM = Laminate1.stressMisesTop
    Ex = np.zeros(Laminate1.nPly)
    E1 = np.zeros(Laminate1.nPly)
    for ii in range(Laminate1.nPly):
        Ex[ii] = np.array([Laminate1.strainLamTopVec[ii][0]])
        E1[ii] = np.array([Laminate1.strainPlyTopVec[ii][0]])
    Laminate1.writeTablePly(
        LaminateCS=True, PlyCS=True, PlyTop=True, PlyBottom=True
    )

    # Initialize lamiante 2 and set parameters
    Laminate2 = laminate()
    Laminate2.LaminateName = 'test1_Laminate2'
    Laminate2.nPly = 11
    Laminate2.theta = [
        -45,
        45,
        90,
        60,
        30,
        0,
        30,
        60,
        90,
        45,
        -45,
    ]
    Laminate2.thkPly = [1] * 11
    Laminate2.E1 = 46000
    Laminate2.E2 = 10000
    Laminate2.G12 = 4600
    Laminate2.nu12 = 0.31
    Laminate2.phi = 0.5
    Laminate2.symmetric = True
    Laminate2.Nx = 100
    Laminate2.Ny = 0
    Laminate2.Nxy = 0
    Laminate2.Mx = 0
    Laminate2.My = 0
    Laminate2.Mxy = 0
    Laminate2.T1 = 1950
    Laminate2.T2 = 48
    Laminate2.C1 = 1480
    Laminate2.C2 = 200
    Laminate2.S12 = 55
    Laminate2.S13 = 55
    Laminate2.S23 = 55
    Laminate2.SymmetricLaminate()
    Laminate2.deg2rad()
    Laminate2.calcPlyPositions()
    Laminate2.calcABD()
    Laminate2.cleanABD()
    Laminate2.calcInverseABD()
    Laminate2.calcLaminateParameters()
    Laminate2.calcStrain()
    Laminate2.calcStress()
    Laminate2.calcMisesStress()
    Laminate2.calcFailureTsaiWu()

    # Plot diagrams and print response tables for laminate 2
    Laminate2.plotStackLayup()
    Laminate2.plotPlyStrainStress()
    Laminate2.writeTablePly(
        LaminateCS=True, PlyCS=False, PlyTop=True, PlyBottom=True
    )

#
