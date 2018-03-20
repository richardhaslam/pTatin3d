Guidelines for Adding Tests
---------------------------

This information relates to the contents within the directory`$PTATIN_DIR/tests`

## Scope

1. Tests should run quickly (ideally, for only a few seconds) and not require more  
resources than those available on a laptop or workstation.

2. While complete coverage is not a current objective, all major features and drivers  
should be covered by at least one test.

## General

1. Each test should be defined in a file `test.py` in a subdirectory of `tests/`.  
See the current tests and [bitbucket.org/dmay/pythontestharness](https://www.bitbucket.org/dmay/pythontestharness)  
for examples and information.

2. Do not include any sub-directories here that do not define tests.

3. Don't forget to add an empty file `__init__.py` within each test sub-directory.

4. Prior to generating any expected output, make sure you execute `make releaseinfo`  
followed by `make all` so that the version will be current in the expected output.

5. For application tests or anything printing residuals, include `-snes_view` so that  
full solver information is available in the expected file.

## Test classes

When defining a new test, please consider which test class / group it should belong to.
Four general classes exist, `unit`, `integration`, `verification`, and `application`.
New tests should be located in the appropriate sub-directory as defined below.

#### a) Unit tests

The sub-directory `tests/unit` contains tests of individual components of the
code, not easily decomposed into smaller units. A good rule of thumb
is that if one of these tests fails, it's obvious where to look for the problem.

#### b) Integration tests

The sub-directory `tests/integration` contains tests which test combinations  
of small numbers of units involved in important tasks. For example, tests of  
linear solves combine  
testing of operators, solvers, etc.

#### c) Verification tests

The sub-directory `tests/verification` contains tests which run benchmarks useful
for verifying the code.

### d) Application tests

The sub-directory `tests/application` contains tests of full application ("driver")  
usage. As opposed to verification tests, these are intended to test the usage  
of the code under conditions which are close to "real" or "production" within  
the constraints on runtime and resource usage mentioned above.

## Executing / Verifying Tests

The testing framework requires python; `python` from your environment is used by default.

### Execute

An individual test can be executed via the following command.

    ./runTests.py -t <class>.<dir1>.<dir2>.<new_test>

e.g.

    ./runTests.py -t unit.json1

The string `<class>.<dir1>.<dir2>.new_test>` defines the complete name of the test.  
To automate and simplify the generation of new tests, the directory structure  
underneath `tests/` is used to define a unique name for each test.

* In the above, `<class>` is either `unit`, `integration`, `verification`, or  
`application`.

* `<dir1>.<dir2>.new_test>` indicates that your test definition can be found in  
`/tests/<class>/dir1/dir2/new_test/test.py`.

For completeness, consider the test defined in  
`$PTATIN_DIR/integration/operator/avx_cuda/compare_mf_3x7x11_p7/test.py`  
To execute this test, one can run

    ./runTests.py -t integration.operator.avx_cuda.compare_mf_3x7x11_p7

### Verify

To verify the output of a test, use the following command

    ./runTests.py -v -t <class>.<dir1>.<dir2>.<new_test>

### Multiple tests

Multiple tests can be executed by providing a comma separated list to the argument `-t`, e.g.

    ./runTests.py -t test1,test2

Similarly, multiple tests can be verified by using `-v` together with comma-separated
list to the argument `-t`, e.g.

    ./runTests.py -v -t test1,test2

### Examine output

Test run in "sandbox" directories named with the test name with `_sandbox`
appended. For example, `unit.json1` will be run from the directory  
`$PTATIN_DIR/tests/unit.json1_sandbox`.  
This directory can be examined for more information on the test run.

Output may be purged with `./runTests.py -p`, or output from a specific test  
may be purged by using the `-t` flag, for example

    ./runTests.py -p unit.json1   # remove unit.json1_sandbox

This purging is automatically performed when a test is (re)run.
