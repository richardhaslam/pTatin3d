Guidelines for adding tests
---------------------------

Do not include any subdirectories here that do not define tests.

Don't forget to add __init__.py in subdirectories here.

Make sure you run `make releaseinfo` and rebuild so that the version will be current
in the expected output.

For application tests or anything printing residuals, include -snes_view so that
full solver information is available in the expected file.
