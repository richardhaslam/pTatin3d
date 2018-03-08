import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 4
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','ptatin_driver_nonlinear_ts.app') + ' -options_file ' + os.path.join(PTATIN_DIR,'src','models','slab_detachment3d','examples','sd3d-pds-test.opts') + ' -options_file ' + os.path.join(thisDir,'more.opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    tol = 1e-6
    t.compareFloatingPoint(re.escape('0 KSP Residual norm'),tol)
    t.compareFloatingPoint(re.escape('1 KSP Residual norm'),tol)
    t.compareFloatingPoint(re.escape('slab base range: gmin'),tol)
    t.compareFloatingPoint(re.escape('slab base range: gmax'),tol)
    t.compareFloatingPoint(re.escape('slab edge range: gmax(x)'),tol)
    t.compareFloatingPoint(re.escape('[step 2] slab edge range: W,Dn'),tol)
    t.compareFloatingPoint(re.escape('slab edge (back face) range: gmax(x)'),tol)
    t.compareFloatingPoint(re.escape('[step 2] slab edge (back face) range: W,Dn'),tol)
    t.compareFloatingPoint(re.escape('slab edge (front face) range: gmin(z)'),tol)
    t.compareFloatingPoint(re.escape('[step 2] slab edge (front face) range: W,Dn'),tol)
    t.compareFloatingPoint(re.escape('topo range:'),tol)
    t.compareFloatingPoint(re.escape('vol,int v.v,vrms:'),tol)

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setVerifyMethod(comparefunc)

  return(t)
