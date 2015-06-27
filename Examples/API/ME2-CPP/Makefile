# If Sherpa-config is not in your $PATH, change this variable to point to your
# Sherpa installation directory manually
SHERPA_PREFIX ?= $(shell Sherpa-config --prefix)

# compiler (eg. w/o MPI c++, w/ MPI mpic++)
CXX    = c++ 
CXXCPP = c++

# add SHERPA_LIBS to LD_LIBRARY PATH environment variable before executing the binary

SHERPA_INCLUDE = $(SHERPA_PREFIX)/include/SHERPA-MC
SHERPA_LIBS    = $(SHERPA_PREFIX)/lib/SHERPA-MC

all: test

test: test.o 
	$(CXX)  test.o -L$(SHERPA_LIBS) -lMEProcess -lToolsOrg -lToolsMath -lToolsPhys -lSherpaMain -lPhasicProcess -o test

test.o: test.cpp
	$(CXXCPP) -fPIC -I$(SHERPA_INCLUDE) -o test.o -c test.cpp

install: test

clean:
	rm -f test test.o
