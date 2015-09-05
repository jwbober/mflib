all: testprime testsqfree testdisc testcrt

CXXFLAGS = -std=c++11 -O2

testprime: testprime.cc ../../include/slint.h
testsqfree: testsqfree.cc ../../include/slint.h
testdisc: testdisc.cc ../../include/slint.h
testcrt: testcrt.cc ../../include/slint.h

clean:
	-rm testprime
	-rm testsqfree
	-rm testdisc
	-rm testcrt
