from LamPy import laminate
import numpy as np

Laminate1 = laminate(nPly=2, theta=np.array([-45.0, 45.0]), thkPly=[1.0, 1.0],
                     E1=6000.0, E2=1000.0, G12=3600.0, nu12=0.32)
Laminate1.LaminateName = "test1_Laminate1"
Laminate1.deg2rad()
Laminate1.calcPlyPositions()
Laminate1.calcABD()
Laminate1.cleanABD()
Laminate1.calcInverseABD()
Laminate1.calcLaminateParameters()
Laminate1.plotStackLayup()
print("Quasi-isotropic laminate +/-45 layup:")
print("Composite laminate Young's modulus in x")
print(Laminate1.ELx)
print("Composite laminate Young's modulus in y")
print(Laminate1.ELy)
