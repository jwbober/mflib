CC = g++
CXX = g++
ROOT = ../..
INCLUDE = $(ROOT)/include
CXXFLAGS = -O2 -std=c++11 -Wall -I/home/bober/include/flint
LDLIBS = -lfftw3 -lflint -larb

OBJECTS = character-info.o \
	  chimaxsum.o \
	  arb_test.o

EXECUTABLES = character-info \
	      chimaxsum \
              arb_test

all: $(EXECUTABLES)

BIN=$(HOME)/bin
links: $(EXECUTABLES)
	ln -f -s $(realpath $(EXECUTABLES)) $(BIN)/

$(OBJECTS): $(INCLUDE)/characters.h $(INCLUDE)/slint.h

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
