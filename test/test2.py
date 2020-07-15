import CLamPy
import numpy as np

nAngle = 101
angles = np.linspace(-90, 90, nAngle)
StressM_Study = [[]]*nAngle
FailureTW_Study = [[]]*nAngle
Strain_xx_Study = [[]]*nAngle
Strain_Fiber_Study = [[]]*nAngle
Laminate1 = CLamPy.laminate()
Laminate1.nPly = 1
Laminate1.thkPly = [2]
for i in range(nAngle):
    Laminate1.theta = [angles[i]]
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
    Strain_xx_Study[i] = Laminate1.strainVec[0][0]
    Strain_Fiber_Study[i] = max(Laminate1.strainPlyTopVec[0][0],
                                Laminate1.strainPlyBotVec[0][0])
    StressM_Study[i] = max(np.max(Laminate1.stressMisesTop),
                           np.max(Laminate1.stressMisesBot))
    FailureTW_Study[i] = np.min(Laminate1.ReserveTsaiWu)
    print("Angle = "+str(round(angles[i], 4))+" deg")
    print("  Failure: "+Laminate1.FailureModeTsaiWu[0])
CLamPy.plotParameterStudyAngle(angles, FailureTW_Study, StressM_Study,
                               Strain_xx_Study, Strain_Fiber_Study, Grid=False,
                               Name="test2")
