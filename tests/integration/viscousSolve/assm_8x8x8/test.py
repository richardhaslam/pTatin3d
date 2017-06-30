import os
import pyTestHarness.unittest as pth
import pyTestHarness.launch as launch
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.split(thisDir)[1]
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_stokes_operators.app') + ' -options_file ' + os.path.join(thisDir,'opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    tol = 1e-6
    t.compareFloatingPoint(re.escape('0 KSP Residual norm'                 ),tol)
    t.compareFloatingPoint(re.escape('1 KSP Residual norm'                 ),tol)
    t.compareFloatingPoint(re.escape('10 KSP Residual norm'                ),tol)
    t.compareFloatingPoint(re.escape('29 KSP Residual norm'                ),tol)

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setVerifyMethod(comparefunc)

  return(t)
