Guidelines for adding tests
---------------------------

This information relates to the contents within the directory`$PTATIN_DIR/tests`

## General 
1. Do not include any sub-directories in here that do not define tests.

2. Don't forget to add an empty file `__init__.py` within each test sub-directories.

3. Prior to generating any expected output, make sure you execute `make releaseinfo` followed by `make all` so that the version will be current
in the expected output.

4. For application tests or anything printing residuals, include -snes_view so that
full solver information is available in the expected file.

## Test classes

When defining a new test, please consider which test class / group it should belong to. Three general classes exist, `unit`, `integration` and `verification`. New tests should be located in the appropriate sub-directory as defined below.

#### a) Unit tests

The sub-directory `tests/unit` contains tests which ...

#### b) Integration tests

The sub-directory `tests/integration` contains tests which ...

#### c) Verification tests
 
The sub-directory `tests/verification` contains tests which ...


## Executing / verifying tests

### Execute

An individual test can be executed via the following command

```python ./runTests.py -t <class>.<dir1>.<dir2>.new_test>```

The string `<class>.<dir1>.<dir2>.new_test>` defines the complete name of the test. Its construction needs further explaination. To automate and simplfy the generation of new tests, the directory struture underneath `tests/` is used to define a unique name for each test.

* In the above, `<class>` is either `unit`, `verification`, or `integration`. 

* `<dir1>.<dir2>.new_test>` indicates that the path to your test definition can be found within `/tests/<class>/dir1/dir2/new_test`.

For completeness, consider the test located here
`integration/operator/avx_cuda/compare_mf_3x7x11_p7`  

To execute this test, we would do the following: 
`python runTests.py -t integration.operator.avx_cuda.compare_mf_3x7x11_p7`


### Verify

To verify the output of a test, use the following command 

```python ./runTests.py -v -t <class.dir1.dir2.new_test>```

### Multiple tests

Multiple tests can be executed by providing a comma separated list to the arguement `-t`, e.g.

```-t test1,test2```

Simarily, multiple tests can be verified by using `-v` together with comma separated list to the arguement `-t`, e.g.

```-v -t test1,test2```