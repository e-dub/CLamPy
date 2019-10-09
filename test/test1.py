from CLamPy import laminate
import numpy as np

plotStackLayup = True
plotPlyStrainStress = True
writeTablePly = True
printTsaiWuData = True
Laminate1 = laminate(nPly=2, theta=np.array([-30, -15]),
                     thkPly=[1, 1.5], E1=38000, E2=9000, G12=3600,
                     nu12=0.32, phi=0.5, symmetric=True,
                     Nx=100, Ny=0, Nxy=0, Mx=0, My=0, Mxy=0, T1=1950, T2=48,
                     C1=1480, C2=200, S12=55, S13=55, S23=55)
Laminate1.LaminateName = "test1_Laminate1"
Laminate1.SymmetricLaminate()
Laminate1.deg2rad()
Laminate1.calcPlyPositions()
Laminate1.calcABD()
Laminate1.cleanABD()
Laminate1.calcInverseABD()
Laminate1.calcLaminateParameters()
Laminate1.calcStrain()
Laminate1.calcStress()
Laminate1.calcMisesStress()
Laminate1.calcFailureTsaiWu()
if plotStackLayup:
    Laminate1.plotStackLayup()
if plotPlyStrainStress:
    Laminate1.plotPlyStrainStress()
if writeTablePly:
    Laminate1.writeTablePly(LaminateCooSys=False, PlyCooSys=True,
                      PlyTop=True, PlyBottom=True)
Laminate2 = laminate(nPly=7, theta=np.array([-45, 45, 0, 30, 60, 90, 45]),
                     thkPly=[1, 1, 2, 2, 5, 2, 3],
                     E1=46000, E2=10000, G12=4600, nu12=0.31,
                     phi=0.5, symmetric=True,
                     Nx=100, Ny=0, Nxy=0, Mx=0, My=0, Mxy=0,
                     T1=1950, T2=48, C1=1480, C2=200, S12=55, S13=55, S23=55)
Laminate2.LaminateName = "test1_Laminate2"
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
if plotStackLayup:
    Laminate2.plotStackLayup()
if plotPlyStrainStress:
    Laminate2.plotPlyStrainStress()
if writeTablePly:
    Laminate2.writeTablePly(LaminateCooSys=True, PlyCooSys=False,
                      PlyTop=True, PlyBottom=True)
if printTsaiWuData:
    print("Laminate 1")
    print("-"*50)
    print("Reserve factor after Tsai-Wu:")
    for i in range(Laminate1.nPly):
        print(Laminate1.ReserveTsaiWu[i])
    print("Failure mode after Tsai-Wu:")
    for i in range(Laminate1.nPly):
        print(Laminate1.FailureModeTsaiWu[i])
    print()
    print("Laminate 2")
    print("-"*50)
    print("Reserve factor after Tsai-Wu:")
    for i in range(Laminate2.nPly):
        print(Laminate2.ReserveTsaiWu[i])
    print("Failure mode after Tsai-Wu:")
    for i in range(Laminate2.nPly):
        print(Laminate2.FailureModeTsaiWu[i])
