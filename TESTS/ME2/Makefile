SOURCES=$(wildcard *.cc)
OBJECTS=$(SOURCES:.cc=.o)
DEPS=$(SOURCES:.cc=.d)
BINS=run_mg

all: $(BINS)

$(BINS): $(OBJECTS)
	$(CXX) -g -o $(BINS) $(OBJECTS)

clean:
	rm -f $(OBJECTS) $(BINS)
