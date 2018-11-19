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
    tol = 1e-13
    t.compareFloatingPoint(re.escape('min[A11_mfo.x-A11_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('max[A11_mfo.x-A11_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('min[A12_mfo.x-A12_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('max[A12_mfo.x-A12_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('min[A21_mfo.x-A21_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('max[A21_mfo.x-A21_asm.x]  = '        ),tol)
    t.compareFloatingPoint(re.escape('min[A_mfo.x-A_asm.x]  = '            ),tol)
    t.compareFloatingPoint(re.escape('max[A_mfo.x-A_asm.x]  = '            ),tol)
    t.compareFloatingPoint(re.escape('min[diagA11_mfo-diagA11_asm]  = '    ),tol)
    t.compareFloatingPoint(re.escape('max[diagA11_mfo-diagA11_asm]  = '    ),tol)
    t.compareFloatingPoint(re.escape('y.y'                                 ),tol)
    t.compareFloatingPoint(re.escape('y2.y2'                               ),tol)

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
