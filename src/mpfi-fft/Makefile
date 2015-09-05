all: test-mpfi-fft mpfi-fft.o

ROOT=/files/math/projects/ntlib
LIB=${ROOT}/lib
INCLUDE=${ROOT}/include

test-mpfi-fft: mpfi-fft.o test-mpfi-fft.o
	gcc -DUSE_MPFI  -o test-mpfi-fft mpfi-fft.o test-mpfi-fft.o -lmpfi -lmpfr
mpfi-fft.o: $(INCLUDE)/mpfi_fft.h mpfi_fft.c
	gcc -c -o mpfi-fft.o mpfi_fft.c -O2
test-mpfi-fft.o: test-mpfi-fft.c $(INCLUDE)/mpfi_fft.h
	gcc -c -o test-mpfi-fft.o test-mpfi-fft.c -O2

clean:
	-rm test-mpfi-fft
	-rm mpfi-fft.o
	-rm test-mpfi-fft.o
