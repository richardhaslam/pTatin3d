import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH') if os.getenv('PETSC_ARCH') else ''
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_stokes_operators.app') + ' -options_file ' + os.path.join(thisDir,'opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    t.compareInteger(re.escape('MatMultA11(MF): iterations'),0) # simply check that 2 iterations were performed

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setVerifyMethod(comparefunc)

  return(t)
