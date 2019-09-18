# If Sherpa-config is not in your $PATH, change this variable to point to your
# Sherpa installation directory manually
SHERPA_PREFIX="/Users/Nappo/work/hep/code/sudakov/sudakov"
# compiler (eg. w/o MPI c++, w/ MPI mpic++)
CXX    = clang++
CXXCPP = clang++

# add SHERPA_LIBS to LD_LIBRARY PATH environment variable before executing the binary

SHERPA_INCLUDE = $(SHERPA_PREFIX)/include/SHERPA-MC
SHERPA_LIBS    = $(SHERPA_PREFIX)/lib/SHERPA-MC

all: check_me

check_me: check_me.o 
	$(CXX)  check_me.o -L$(SHERPA_LIBS) -lMEProcess -lToolsOrg -lToolsMath -lToolsPhys -lSherpaMain -lPhasicProcess -o check_me

check_me.o: check_me.cpp
	$(CXXCPP) -std=c++11 -stdlib=libc++ -fPIC -I$(SHERPA_INCLUDE) -o check_me.o -c check_me.cpp

install: check_me

clean:
	rm -f check_me check_me.o
