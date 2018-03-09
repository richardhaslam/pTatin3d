#!/usr/bin/env python
################################################################################
# Helper script to generate test subsets                                       #
#                                                                              #
# Prints a comma-separated list of tests, to be used with pyTestHarness's      #
# -t option.                                                                   #
#                                                                              #
# Run with -h for information on options.                                      #
#                                                                              #
# Looks for files called 'test.py' and                                         #
# assumes that they define tests named by the relative path to their enclosing #
# directory, from $PTATIN_DIR, with directory separators replaced with dots.   #
################################################################################

import os
import sys
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--smoke',action='store_true',help='Override other arguments and run smoke tests')
parser.add_argument('--skip',help='comma-separated list of subdirectories to skip, e.g. --skip cuda,opencl')
parser.add_argument('--prefix',help='subtree of PTATIN_DIR to restrict to, e.g. -- prefix unit/json')
args = parser.parse_args()

# Custom logic for smoke test
if (args.smoke) :
    print('integration.operator.compare_mf_asm_8x9x11,integration.viscousSolve.asm_8x8x8,application.linear.sinker_galerkin_8x8x8')
    sys.exit(0)

# Look for subdirectories of $PTATIN_DIR/tests
PTATIN_DIR = os.getenv('PTATIN_DIR')
if not PTATIN_DIR :
    print("You must define PTATIN_DIR in your environment. Exiting.")
    sys.exit(1)

# Directory names to ignore
ignore = set(args.skip.split(','))
if '' in ignore :
    ignore.remove('')

# Subtree to search in
if args.prefix :
    prefix = args.prefix
else :
    prefix = ''

# Walk the (sub)tree to create a comma-separated list of tests
resultString=  ''
testsDir = os.path.join(PTATIN_DIR,'tests')
for (root, dirs, files) in os.walk(os.path.join(testsDir,prefix)) :
    if 'test.py' in files :
        testName = os.path.relpath(root,testsDir).replace(os.sep,'.')
        if set(testName.split('.')).isdisjoint(ignore) :
            if (len(resultString) > 0) :
                resultString += ','
            resultString += testName
print(resultString) 
