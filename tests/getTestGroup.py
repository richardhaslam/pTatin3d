#!/usr/bin/env python
################################################################################
# Helper script to generate test subsets                                       #
#                                                                              #
# Prints a comma-separated list of tests, to be used with pyTestHarness's      #
# -t option.                                                                   #
#                                                                              #
# The first argument can be either a custom group name (e.g. smoke), or a      #
# subdirectory of $PTATIN_DIR/tests. It is used to choose a subset of tests.   #
# This may also be the special value "all"                                     #
#                                                                              #
# Additional arguments can be supplied of the form "no_xxx", which will ignore #
# any directory with name "xxx". This is useful to skip things like "cuda"     #
# directories when the library is not configured with CUDA                     #
#                                                                              #
# When a subdirectory is specified, looks for files called 'test.py' and       #
# assumes that they define tests named by the relative path to their enclosing #
# directory, from $PTATIN_DIR.                                                 #
################################################################################

import os
import sys

# Custom Groups
if (len(sys.argv) > 1) :
    if 'smoke' == sys.argv[1] :
        print('integration.operator.compare_mf_asm_8x9x11,integration.viscousSolve.asm_8x8x8,application.linear.sinker_galerkin_8x8x8')
        if (len(sys.argv) > 2) :
            print('Warning:getTestGroup.py ignoring additional arguments')
        sys.exit(0)
    if 'classic' == sys.argv[1] :
        print('integration.operator.compare_mf_asm_8x9x11,integration.viscousSolve.asm_8x8x8,application.linear.sinker_galerkin_8x8x8,application.linear.sinker_hybrid_8x8x8,application.linear.sinker_hybrid_8x8x8_p8,application.linear.sinker_mf_8x8x8')
        if (len(sys.argv) > 2) :
            print('Warning:getTestGroup.py ignoring additional arguments')
        sys.exit(0)
    if 'all' == sys.argv[1] :
        prefix = ''
    else :
        prefix = sys.argv[1]
else :
    prefix = ''

# If no custom group matched, look for subdirectories of $PTATIN_DIR/tests
PTATIN_DIR = os.getenv('PTATIN_DIR')
if not PTATIN_DIR :
    print("You must define PTATIN_DIR in your environment. Exiting.")
    sys.exit(1)

resultString=  ''
testsDir = os.path.join(PTATIN_DIR,'tests')
ignore = set([entry[3:] for entry in sys.argv[2:]])            # strip off "no_"
for (root, dirs, files) in os.walk(os.path.join(testsDir,prefix)) :
    if 'test.py' in files :
        testName = os.path.relpath(root,testsDir).replace(os.sep,'.')
        if set(testName.split('.')).isdisjoint(ignore) :
            if (len(resultString) > 0) :
                resultString += ','
            resultString += testName
print(resultString) 
