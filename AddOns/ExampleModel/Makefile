# If Sherpa-config is not in your $PATH, change this variable to point to your
# Sherpa installation directory manually
SHERPA_PREFIX=$(Sherpa-config --prefix)


SOURCES=SM_Zprime.C Interaction_Model_SM_Zprime.C

CC=g++
CXXFLAGS=-g -O2 -pedantic -Wall -fPIC -DPIC

CXXFLAGS+= -I$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --incdir)
LDFLAGS=$(shell $(SHERPA_PREFIX)/bin/Sherpa-config --ldflags)

OBJECTS=$(SOURCES:.C=.o)

all: $(SOURCES) $(OBJECTS) 
	$(CC) -shared -o libSMZprime.so $(OBJECTS) $(LDFLAGS)

.C.o:
	$(CC) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) libSMZprime.so
