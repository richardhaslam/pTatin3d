import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH') if os.getenv('PETSC_ARCH') else ''
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_linear_ts.app') + ' -options_file ' + os.path.join(PTATIN_DIR,'src','models','analytics_vv','examples','analytics-solcx-isoviscous.opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    t.compareFloatingPoint(re.escape('1 KSP Component U,V,W,P residual norm '),1e-6)
    t.compareFloatingPoint(re.escape('[inorms] '),1e-6)

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
