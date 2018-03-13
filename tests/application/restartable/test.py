import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH') if os.getenv('PETSC_ARCH') else ''
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 4
  launch = [\
          os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_ts_init.app') + ' -options_file ' + os.path.join(thisDir,'opts') + ' -options_file ' + os.path.join(thisDir,'more1.opts'),\
          os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_ts_init.app') + ' -options_file ' + os.path.join(thisDir,'opts') + ' -options_file ' + os.path.join(thisDir,'more2.opts'),\
          os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_ts_init.app') + ' -options_file ' + os.path.join(thisDir,'opts') + ' -options_file ' + os.path.join(thisDir,'more3.opts'),\
          ]
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    tol = 1e-10
    # Note - These KSP output tests are not ideal. These will only look at the U component, and tol is absolute
    t.compareFloatingPoint(re.escape('0 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('1 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('2 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('3 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('4 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('5 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('6 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('7 KSP Component U,V,W,P residual norm ['),tol)
    t.compareFloatingPoint(re.escape('8 KSP Component U,V,W,P residual norm ['),tol)

    t.compareFloatingPoint(re.escape('0 SNES Function norm'),tol)
    t.compareFloatingPoint(re.escape('1 SNES Function norm'),tol)

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
