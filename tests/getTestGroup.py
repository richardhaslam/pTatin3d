#!/usr/bin/env python
################################################################################
# Helper script to generate test subsets                                       #
#                                                                              #
# Supply a single argument, which is either a custom group name (e.g. smoke),  #
# or is a subdirectory of $PTATIN_DIR/tests.                                   #
#                                                                              #
# Prints a comma-separated list of tests, to be used with pyTestHarness's      #
# -t option.                                                                   #
#                                                                              #
# When a subdirectory is specified, looks for files called 'test.py' and       #
# assumes that they define tests named by the relative path to their enclosing #
# directory, from $PTATIN_DIR.                                                 #
################################################################################

import os
import sys

# Custom Groups
if 'smoke' in sys.argv :
    print('integration.operator.compare_mf_asm_8x9x11,integration.viscousSolve.asm_8x8x8,application.linear.sinker_galerkin_8x8x8')
    sys.exit(0)

if 'classic' in sys.argv :
    print('integration.operator.compare_mf_asm_8x9x11,integration.viscousSolve.asm_8x8x8,application.linear.sinker_galerkin_8x8x8,application.linear.sinker_hybrid_8x8x8,application.linear.sinker_hybrid_8x8x8_p8,application.linear.sinker_mf_8x8x8')
    sys.exit(0)

# If no custom group matched, look for subdirectories of $PTATIN_DIR/tests
PTATIN_DIR = os.getenv('PTATIN_DIR')
if not PTATIN_DIR :
    print("You must define PTATIN_DIR in your environment. Exiting.")
    sys.exit(1)

prefix = sys.argv[1]
resultString=  ""
testsDir = os.path.join(PTATIN_DIR,'tests')
for (root, dirs, files) in os.walk(os.path.join(testsDir,prefix)) :
    if 'test.py' in files :
        testName = os.path.relpath(root,testsDir).replace(os.sep,'.')
        if (len(resultString) > 0) :
            resultString += ","
        resultString += testName
print(resultString) 
