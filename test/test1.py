from CLamPy import laminate
import numpy as np

# Initialize lamiante 1 and set parameters
Laminate1 = laminate()
Laminate1.LaminateName = "test1_Laminate1"
Laminate1.nPly = 2
Laminate1.theta = np.array([-30, -15])
Laminate1.thkPly = [1, 1.5]
Laminate1.E1 = 38000
Laminate1.E2 = 9000
Laminate1.G12 = 3600
Laminate1.nu12 = 0.32
Laminate1.phi = 0.5
Laminate1.symmetric = True
Laminate1.Nx = 100
Laminate1.Ny = 0
Laminate1.Nxy = 0
Laminate1.Mx = 0
Laminate1.My = 0
Laminate1.Mxy = 0
Laminate1.T1 = 1950
Laminate1.T2 = 48
Laminate1.C1 = 1480
Laminate1.C2 = 200
Laminate1.S12 = 55
Laminate1.S13 = 55
Laminate1.S23 = 55

# Calculate values for laminate 1
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

# Plot laminate 1
Laminate1.plotStackLayup()
Laminate1.plotPlyStrainStress()
Laminate1.writeTablePly(LaminateCS=False, PlyCS=True, PlyTop=True,
                        PlyBottom=True)

# Initialize lamiante 2 and set parameters
Laminate2 = laminate()
Laminate2.LaminateName = "test1_Laminate2"
Laminate2.nPly = 7
Laminate2.theta = [-45, 45, 0, 30, 60, 90, 45]
Laminate2.thkPly = [1, 1, 2, 2, 5, 2, 3]
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
Laminate2.writeTablePly(LaminateCS=True, PlyCS=False, PlyTop=True,
                        PlyBottom=True)

# Print Tsai-Wu data for laminates 1 and 2
print()
print("Laminate 1")
print("-"*50)
print("Reserve factor after Tsai-Wu:")
for i in range(Laminate1.nPly):
    print("  Ply "+str(i)+": "+str(round(Laminate1.ReserveTsaiWu[i], 4)))
print("Failure mode after Tsai-Wu:")
for i in range(Laminate1.nPly):
    print("  Ply "+str(i)+": "+Laminate1.FailureModeTsaiWu[i])
print()
print("Laminate 2")
print("-"*50)
print("Reserve factor after Tsai-Wu:")
for i in range(Laminate2.nPly):
    print("  Ply "+str(i)+": "+str(round(Laminate2.ReserveTsaiWu[i], 4)))
print("Failure mode after Tsai-Wu:")
for i in range(Laminate2.nPly):
    print("  Ply "+str(i)+": "+Laminate2.FailureModeTsaiWu[i])
