CC = g++
CXX = g++
ROOT = ../..
INCLUDE = $(ROOT)/include
CXXFLAGS = -O2 -std=c++11 -g -I/home/bober/include/flint -isystem /usr/include/eigen3 -fpermissive -Wall -Wno-sign-compare -DUSE_ARB -Wno-unused-function
CFLAGS = -O2 -g -I/home/bober/include/flint -Wall
ifdef STATIC_LIB_DIR
    LDLIBS = -lpthread -lfftw3 $(STATIC_LIB_DIR)/libarb.a $(STATIC_LIB_DIR)/libflint.a $(STATIC_LIB_DIR)/libntl.a $(STATIC_LIB_DIR)/libmpfr.a $(STATIC_LIB_DIR)/libgmp.a
else
    LDLIBS = -L$(HOME)/lib -lfftw3 -lflint -larb -lgmp -lntl -lpthread
endif


OBJECTS = modform_cc.o \
	  modform_modp.o \
	  modp_newbasis.o \
	  modp_trace.o \
	  classnumbers.o \
	  newforms.o \
	  print-traces.o \
	  modp_newforms.o \
	  modform_acb.o \
	  newforms_acb.o \
	  arb_extras.o \
	  identify-fields.o \
	  hecke-polynomial.o \
	  modp_basis.o \
	  weight1.o

EXECUTABLES = modp_newbasis \
	      modp_trace \
	      newforms \
	      print-traces \
	      modp_newforms \
	      newforms_acb \
	      identify-fields \
	      hecke-polynomial \
	      modp_basis \
	      weight1 \
	      convert-to-binary \
	      cat-binary-format

all: $(EXECUTABLES)

BIN=$(HOME)/bin
links: $(EXECUTABLES)
	ln -f -s $(realpath $(EXECUTABLES)) $(BIN)/

$(OBJECTS): $(INCLUDE)/characters.h $(INCLUDE)/slint.h modform_modp.h modform_cc.h $(INCLUDE)/characters.h
arb_extras.o: arb_extras.c
	gcc -c $(CFLAGS) arb_extras.c -o arb_extras.o
classnumbers.o: classnumbers.cc
identify-fields.o: identify-fields.cc $(INCLUDE)/slint.h $(INCLUDE)/characters.h

modp_newbasis: modform_modp.o modp_newbasis.o classnumbers.o
modp_basis: modform_modp.o modp_basis.o classnumbers.o
modp_trace: modform_modp.o modp_trace.o classnumbers.o
weight1: modform_modp.o weight1.o classnumbers.o
newforms: newforms.o modform_cc.o modform_modp.o classnumbers.o
newforms_acb: arb_extras.o newforms_acb.o modform_acb.o modform_modp.o classnumbers.o modform_cc.o
modp_newforms: modp_newforms.o modform_modp.o classnumbers.o
print-traces: print-traces.o modform_cc.o modform_modp.o classnumbers.o
identify-fields: identify-fields.o
hecke-polynomial: modform_modp.o hecke-polynomial.o classnumbers.o
modform_acb.o: modform_acb.h
cat-binary-format.o: mfformat.h
convert-to-binary: mfformat.h

clean:
	-rm $(OBJECTS)
	-rm $(EXECUTABLES)
