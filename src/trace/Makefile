CC = g++
CXX = g++
ROOT = ../..
INCLUDE = $(ROOT)/include
CXXFLAGS = -O2 -std=c++11 -g -I/home/bober/include/flint -fpermissive -Wall
LDLIBS = -lfftw3 -lflint

OBJECTS = trace-formula.o \
	  weight1.o \
	  test-tf.o \
	  trace-unsieved.o \
	  trace-sieved.o \
	  trace-TnTm.o \

EXECUTABLES = weight1 \
	      test-tf \
	      trace-unsieved \
	      trace-sieved \
	      trace-TnTm \

all: $(EXECUTABLES)

BIN=$(HOME)/bin
links: $(EXECUTABLES)
	ln -f -s $(realpath $(EXECUTABLES)) $(BIN)/

$(OBJECTS): $(INCLUDE)/characters.h $(INCLUDE)/slint.h trace-formula.h

write-classnumbers.o: write-classnumbers.cc
write-classnumbers: write-classnumbers.o
	g++ -o write-classnumbers write-classnumbers.o -lpari

weight1: trace-formula.o weight1.o
test-tf: trace-formula.o test-tf.o
trace-unsieved: trace-formula.o trace-unsieved.o
trace-sieved: trace-formula.o trace-sieved.o
trace-TnTm: trace-formula.o trace-TnTm.o

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
	-rm write-classnumbers
	-rm write-classnumbers.o
