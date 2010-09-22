# Compile via 'make SHERPA_PREFIX=/path/to/sherpa'
# If using the 'heracmn' alone, add option 'HERACMN=1'
# If using the H1 interface, add option 'H1=1'

LIBNAME = libSherpaHEPEVT
SOURCES = HEPEVT_Interface.C

CC = g++
CXXFLAGS = -g -O2 -pedantic -Wall -fPIC -DPIC

ifeq ($(H1),1)
CXXFLAGS += -DUSING__H1_Mode
else
ifeq ($(HERACMN),1)
CXXFLAGS += -DUSING__HERACMN
endif
endif

CXXFLAGS += -I$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --incdir)

OBJECTS=$(SOURCES:.C=.o)

all: $(SOURCES) $(OBJECTS)
	$(CC) -shared -o $(LIBNAME).so $(OBJECTS)

.C.o:
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) $(LIBNAME).so
