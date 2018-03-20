import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH') if os.getenv('PETSC_ARCH') else ''
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_linear_ts.app') + ' -options_file ' + os.path.join(PTATIN_DIR,'src','models','static_box','examples','staticbox-ex1.opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    t.compareLiteral(re.escape('[staticBox] ||v0 - v0_exact||_inf <= 1.0e-12 '))
    t.compareLiteral(re.escape('[staticBox] ||v1 - v1_exact||_inf <= 1.0e-12 '))
    t.compareLiteral(re.escape('[staticBox] ||v2 - v2_exact||_inf <= 1.0e-12 '))
    t.compareLiteral(re.escape('[staticBox] ||p - p_exact||_inf <= 1.0e-12 '  ))
    t.compareLiteral(re.escape('[staticBox] rho = '                           ))

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
