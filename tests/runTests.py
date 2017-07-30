#!/usr/bin/env python
#################################################################################
#                                                                               #
#                                pTatin3D Tests                                 #
#                                                                               #
#################################################################################
# Usage:                                                                        #
#  ./runTests.py [-h] [-other_pyTestHarness_options]                            #
#                                                                               #
# Design principles:                                                            #
# * It should be easy to add, move, or rename tests                             #
# * It should be easy to run a meaningful subset of tests                       #
# * Tests should have descriptive names                                         #
#                                                                               #
#################################################################################

import os
import sys

# PTATIN_DIR and PETSC_ARCH are used to make easily-copyable tests
thisDir = os.path.dirname(os.path.abspath(__file__))
PTATIN_DIR = os.getenv('PTATIN_DIR')
if not PTATIN_DIR :
    print("You must define PTATIN_DIR in your environment. Exiting.")
    sys.exit(1)
if os.path.join(PTATIN_DIR,'tests') != thisDir :
    print("PTATIN_DIR is not set correctly in your environment. Exiting.")
    sys.exit(1)
PETSC_ARCH = os.getenv('PETSC_ARCH')
if not PETSC_ARCH :
    print("You must define PETSC_ARCH in your environment. Exiting.")
    sys.exit(1)

# bitbucket.org/dmay/pythontestharness
sys.path.append(os.path.join(PTATIN_DIR,'tests','pythontestharness','lib'))  # overrides
try :
  import numpy
except Exception as errorMessage :
  print(format(errorMessage))
  print("********************")
  print("numpy is required to use the required pyTestHarness module.")
  print("Make sure that numpy is available with your python installation.")
  print("********************")
  sys.exit(1)
try :
  import pyTestHarness.harness as pyth_harness
except Exception as errorMessage :
  if not sys.exc_info()[-1].tb_next :
    print(format(errorMessage))
    print("********************")
    print("The required python library pyTestHarness was not found. Exiting.")
    print("If pyTestHarness is installed on your system, ensure pythontestharness/lib is included in the environment variable PYTHONPATH.")
    print("If pyTestHarness is not installed, obtain the source by executing the following:")
    print("  git clone https://bitbucket.org/dmay/pythontestharness " + os.path.join(PTATIN_DIR,'tests','pythontestharness'))
    print("********************")
    sys.exit(1)

#  Interpret any child directory containing "test.py" as defining as a test
#  with a name defined as the relative path to the containing directory.
allTests = []
for (root, dirs, files) in os.walk(thisDir) :
    if 'test.py' in files :
        relpath = os.path.relpath(root,thisDir)
        mod = relpath.replace(os.sep,'.') + '.test'
        exec('import ' + mod)
        exec('allTests.append('+mod+'.test())')

# Run tests
os.environ['PYTHONUNBUFFERED'] = str('1')
h = pyth_harness.pthHarness(allTests)
h.execute()
h.verify()
