# Sherpa-config has to be in your $PATH for the following to work
all:
	g++ -I`Sherpa-config --incdir` -fPIC -c Userhook_Example.C -o Userhook_Example.o
	g++ -o libSherpaUserhookExample.so Userhook_Example.o `Sherpa-config --ldflags` -shared
