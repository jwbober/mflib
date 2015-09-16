CC = g++
CXX = g++
ROOT = ../..
INCLUDE = $(ROOT)/include
CXXFLAGS = -O2 -std=c++11 -Wall
LDLIBS = -lfftw3

OBJECTS = character-info.o \

EXECUTABLES = character-info \

all: $(EXECUTABLES)

BIN=$(HOME)/bin
links: $(EXECUTABLES)
	ln -f -s $(realpath $(EXECUTABLES)) $(BIN)/

$(OBJECTS): $(INCLUDE)/characters.h $(INCLUDE)/slint.h

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
