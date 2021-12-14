# list test process names; it is expected that for each process name, a
# corresponding file `name.yaml` exists within this directory, which contains
# the Sherpa setup to run for the test; note that a special file `Base.yaml`
# exists which will also be read in by Sherpa, in order to provide
# configuration options common to all tests
TEST_PROCS = eePP eeWW eeZP eeZZ eebb eedd eemumu eett eeuu duWZ duWP uuPG udWG uuZG

# define VERBOSE, such that the test suite log file is printed when tests fail;
# this is useful in particular to diagnose issues that happen in the CI only
VERBOSE = 1

# define tests and the script used to run each test (via Automake)
TESTS = $(TEST_PROCS:=.yaml)
LOG_COMPILER = bash run.sh

# the run script runs each test in a directory named after the process;
# Automake doesn't know about this output, so add it manually to `make clean`
clean-local:
	-rm -rf $(TEST_PROCS)
