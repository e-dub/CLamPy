import CLamPy
import numpy as np

angles = np.linspace(-90, 90, 101)
StressM_Study = [[]]*len(angles)
FailureTW_Study = [[]]*len(angles)
Strain_xx_Study = [[]]*len(angles)
Strain_Fiber_Study = [[]]*len(angles)
Laminate1 = pyCLT.laminate(nPly=1, thkPly=[2])
for i in range(len(angles)):
    Laminate1.theta=[angles[i]]
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
    print(Laminate1.FailureModeTsaiWu)
pyCLT.plotParameterStudyAngle(angles, FailureTW_Study, StressM_Study,
                              Strain_xx_Study, Strain_Fiber_Study, Grid=False,
                              Name="test2")
