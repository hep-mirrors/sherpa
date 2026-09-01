# Determine appropriate extension for plugin library
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  SHLIB_EXT := dylib
else
  SHLIB_EXT := so
endif

# Sherpa-config has to be in your $PATH for the following to work
LIBNAME := libSherpaUserhookExample.$(SHLIB_EXT)
SRC := Userhook_Example.C
CXX := $(shell Sherpa-config --cxx)
CXXFLAGS := $(shell Sherpa-config --cppflags) -fPIC -Wall
LDFLAGS := -shared $(shell Sherpa-config --ldflags)

all: $(LIBNAME)

$(LIBNAME): $(SRC)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o *.so *.dylib

