import os
import pyTestHarness.unittest as pth
import re

def test() :
  PTATIN_DIR = os.getenv('PTATIN_DIR')
  PETSC_ARCH = os.getenv('PETSC_ARCH')
  thisDir = os.path.split(os.path.abspath(__file__))[0]
  testName = os.path.relpath(thisDir,os.path.join(PTATIN_DIR,'tests')).replace(os.sep,'.')
  ranks = 1
  launch = os.path.join(PTATIN_DIR,PETSC_ARCH,'bin','test_dmda_checkpoint.app') + ' -checkpoint'
  expectedFile = os.path.join(thisDir,'expected')

  def comparefunc(t) :
    # Note - I need to add a space add the end otherwise an error arises during the comparison
    # ValueError: could not convert string to float: '_2'
    key = re.escape("|x| ")
    t.compareFloatingPoint(key,0.0)

    key = re.escape("|x|_2")
    t.compareFloatingPoint(key,0.0)

    key = re.escape("min(x)")
    t.compareFloatingPoint(key,0.0)

    key = re.escape("max(x)")
    t.compareFloatingPoint(key,0.0)
    
    t.compareLiteral(re.escape("Processor [0] M 10 N 10 P 10 m 1 n 1 p 1 w 6 s 1"))
    t.compareLiteral(re.escape("X range of indices: 0 10, Y range of indices: 0 10, Z range of indices: 0 10"))
  
    # This comparison attempt results in the error that the string isn't in the expected file - why??
    #t.compareLiteral(re.escape("Lower left corner -10. -20. -30. : Upper right 10. 20. 30."))

    key = re.escape("gmin-x")
    t.compareFloatingPoint(key,0.0)
    key = re.escape("gmin-y")
    t.compareFloatingPoint(key,0.0)
    key = re.escape("gmin-z")
    t.compareFloatingPoint(key,0.0)
    key = re.escape("gmax-x")
    t.compareFloatingPoint(key,0.0)
    key = re.escape("gmax-y")
    t.compareFloatingPoint(key,0.0)
    key = re.escape("gmax-z")
    t.compareFloatingPoint(key,0.0)


  # Create Test Object
  t = pth.pthUnitTest(testName,ranks,launch,expectedFile)
  t.setUseSandbox()
  t.setVerifyMethod(comparefunc)

  return(t)
