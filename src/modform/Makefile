CC = g++
CXX = g++
ROOT = ../..
INCLUDE = $(ROOT)/include
CXXFLAGS = -O2 -std=c++11 -g -I/home/bober/include/flint -isystem /usr/include/eigen3 -fpermissive -Wall
LDLIBS = -L$(HOME)/lib -lfftw3 -lflint -larb

OBJECTS = modform_cc.o \
	  modform_modp.o \
	  modp_newbasis.o \
	  modp_trace.o \
	  classnumbers.o \
	  newforms.o

EXECUTABLES = modp_newbasis \
	      modp_trace \
	      newforms \

all: $(EXECUTABLES)

BIN=$(HOME)/bin
links: $(EXECUTABLES)
	ln -f -s $(realpath $(EXECUTABLES)) $(BIN)/

$(OBJECTS): $(INCLUDE)/characters.h $(INCLUDE)/slint.h modform_modp.h modform_cc.h

classnumbers.o: classnumbers.cc

modp_newbasis: modform_modp.o modp_newbasis.o classnumbers.o
modp_trace: modform_modp.o modp_trace.o classnumbers.o
newforms: newforms.o modform_cc.o modform_modp.o classnumbers.o

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
