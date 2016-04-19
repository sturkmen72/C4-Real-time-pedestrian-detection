CC=g++
WLEVEL=-Wall
CFLAGS=-c $(WLEVEL) -march=native -O3 `pkg-config --cflags opencv` -ffast-math
LDFLAGS= `pkg-config --libs opencv`
SOURCES= Pedestrian_ICRA.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)

all:  detect

depend:
	g++ $(CFLAGS) -MM $(SOURCES) > .deps

clean:
	rm -rf *.o detect *~ *.bak

detect: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $*.cpp

include .deps
