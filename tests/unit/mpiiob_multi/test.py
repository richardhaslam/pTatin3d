import os
import re
import pyTestHarness.test as pthtest

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 4

  # Launch the executable twice
  launch = [\
          os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_mpiiob.app') + ' -test_id 3 -data_len 1000',\
          os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_mpiiob.app') + ' -test_id 4 -data_len 1000',\
          ]
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    # We expect an exact match
    t.compareLiteral(re.escape('++ [multi]'));

  # Create Test Object
  t = pthtest.Test(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
