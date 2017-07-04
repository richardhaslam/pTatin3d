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

import pyTestHarness.harness as pyth_harness   # bitbucket.org/dmay/pyTestHarness
import os

# PTATIN_DIR and PETSC_ARCH are used to make easily-copyable tests
thisDir = os.path.dirname(os.path.abspath(__file__))
PTATIN_DIR = os.getenv('PTATIN_DIR')
if not PTATIN_DIR :
    print("You must define PTATIN_DIR in your environment. Exiting.")
    exit(1)
if os.path.join(PTATIN_DIR,'tests') != thisDir :
    print("PTATIN_DIR is not set correctly in your environment. Exiting.")
    exit(1)
PETSC_ARCH = os.getenv('PETSC_ARCH')
if not PETSC_ARCH :
    print("You must define PETSC_ARCH in your environment. Exiting.")
    exit(1)

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
# TODO: these need to be somehow sandboxed, so that output files (even if you don't think you're supposed to generate any) don't clash
os.environ['PYTHONUNBUFFERED'] = str('1')
launcher = pyth_harness.pthHarness(allTests)
launcher.execute()
exitcode = launcher.verify()
launcher.clean()
exit(exitcode)
