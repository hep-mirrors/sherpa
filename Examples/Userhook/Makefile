SHERPA_PREFIX ?= $(shell Sherpa-config --prefix)
SOURCES=Userhook_Example.C

CXXFLAGS+= -I$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --incdir) -fPIC
LDFLAGS=$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --ldflags) -shared

OBJECTS=$(SOURCES:.C=.o)

all: $(SOURCES) $(OBJECTS) 
	$(CXX) -o libSherpaUserhookExample.so $(OBJECTS) $(LDFLAGS)

.C.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@
