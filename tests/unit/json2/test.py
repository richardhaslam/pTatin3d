import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_cjson.app') + ' -test_id 2'
  expectedFile = os.path.join(thisDir,'fat_test_geom.json.expected')

  def comparefunc(t) :
    t.compareUnixDiff()

  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setVerifyMethod(comparefunc)
  t.setComparisonFile('fat_test_geom.json')

  return(t)
