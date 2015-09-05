#test-values test-primitive test-conductor test-character-extension

CC=g++
ROOT=/files/math/projects/ntlib
LIB=${ROOT}/lib
INCLUDE=${ROOT}/include
CXXFLAGS = -std=c++11 -O2

LDLIBS = -lfftw3 -lmpfi -lmpfr

CXXFLAGS = -std=c++11 -O2
OBJECTS = dft-test.o test-primitive.o test-values.o test-conductor.o test-character-extension.o
EXECUTABLES = dft-test test-primitive test-values test-conductor test-character-extension

all: $(EXECUTABLES)
$(OBJECTS): $(INCLUDE)/characters.h

dft-test: dft-test.o ${LIB}/mpfi-fft.o
dft-test.o : CXXFLAGS += -DUSE_MPFI
dft-test.o: dft-test.cc ../../include/characters.h

#dft-test: dft-test.cc characters.h ../mpfi-fft/mpfi-fft.o
#	g++ -DUSE_MPFI -O3 -o dft-test dft-test.cc ../mpfi-fft/mpfi-fft.o -lfftw3 -lmpfi -lmpfr -std=c++11
#test-primitive: characters.h test-primitive.cc
#	g++ -O3 -o test-primitive test-primitive.cc -lfftw3 -lmpfr -std=c++11
#test-values: characters.h test-values.cc
#	g++ -O3 -o test-values test-values.cc -lfftw3 -lmpfr -std=c++11
#test-conductor: test-conductor.cc characters.h
#	g++ -std=c++11 -O3 -o test-conductor test-conductor.cc -pthread -lfftw3 -std=c++11
#test-character-extension: test-character-extension.cc characters.h
#	g++ -std=c++11 -O3 -o test-character-extension test-character-extension.cc -pthread -lfftw3 -std=c++11

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
