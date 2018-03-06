import os
import re
import pyTestHarness.unittest as pth

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 4
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_databucket_checkpoint.app') + ' -options_file ' + os.path.join(thisDir,'opts')
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    t.compareLiteral('coord statistic 0 PASSED')
    t.compareLiteral('coord statistic 1 PASSED')
    t.compareLiteral('coord statistic 2 PASSED')
    t.compareLiteral('coord statistic 3 PASSED')
    t.compareLiteral('coord statistic 4 PASSED')
    t.compareLiteral('coord statistic 5 PASSED')
    t.compareLiteral('coord statistic 6 PASSED')
    t.compareLiteral('coord statistic 7 PASSED')
    t.compareLiteral('coord statistic 8 PASSED')

  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
