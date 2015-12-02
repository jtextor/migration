
CFLAGS=-O2
CXXFLAGS=-O2
CXX=g++
SIMHEADERS=jacobi.h params.h
OBJECTS=INIReader.o activation-time.o ini.o

all: activation-time residence-times

activation-time: $(OBJECTS)
	g++ $^ -o $@

residence-times: residence-times.cpp $(SIMHEADERS)
	g++ $(CXXFLAGS) $< -o $@ 

clean:
	rm -f $(OBJECTS) activation-times residence-times
