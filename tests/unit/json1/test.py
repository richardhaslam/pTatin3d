import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH') if os.getenv('PETSC_ARCH') else ''
  PETSC_ARCH = PETSC_ARCH if PETSC_ARCH else ''
  print(PETSC_ARCH)
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_cjson.app') + ' -test_id 1 -filename1 ' + os.path.join(PTATIN_DIR,'src','tests','test_geom.json')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    t.compareUnixDiff()

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  # Ignore ** (header) lines
  t.appendKeywords('**')

  return(t)
